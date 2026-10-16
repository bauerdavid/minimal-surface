"""Add the licenses of the libraries bundled into a wheel to the wheel itself.

The wheels carry ITK (and what it depends on), the OpenMP runtime and SimpleITK, and the extension
is compiled against Eigen. Their licenses ask for their notices to travel with the binaries, so
this collects them into THIRD_PARTY_LICENSES.txt and stores that file in the wheel's
<name>.dist-info/licenses/ directory, where PyPI and `pip` metadata tools expect license files.

    licenses.py conda --prefix ENV --wheel W --eigen DIR   Windows and macOS: the libraries came from
                                                           conda packages, whose license files are used
    licenses.py linux --wheel W --sdk-licenses DIR --eigen DIR
                                                           Linux: ITK and SimpleITK are statically linked
                                                           (licenses gathered when the SDK was built) and
                                                           the OpenMP runtime comes from the system

Both fail when the wheel holds a library they cannot account for, so a new dependency cannot slip
into a release without its license.
"""
import argparse
import base64
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile

NAME = "THIRD_PARTY_LICENSES.txt"
LIBRARY_SUFFIXES = (".dll", ".dylib", ".so")

# Distributed with Visual Studio rather than as conda packages.
MICROSOFT_RUNTIME = re.compile(r"^(vcomp|msvcp|vcruntime|concrt|vccorlib)\d*(_\d+)?\.dll$", re.I)
MICROSOFT_NOTE = (
    "vcomp140.dll (the OpenMP runtime), msvcp140.dll and the other Visual C++ runtime DLLs in this\n"
    "wheel are part of the Microsoft Visual C++ Redistributable, redistributed under Microsoft's terms\n"
    "for distributable code in Visual Studio. See\n"
    "https://learn.microsoft.com/en-us/cpp/windows/redistributing-visual-cpp-files\n"
)


def bundled_libraries(wheel):
    """The names of the shared libraries copied into the wheel, with the hash the repair tool added removed."""
    names = []
    with zipfile.ZipFile(wheel) as z:
        for member in z.namelist():
            base = os.path.basename(member)
            if not base.lower().endswith(LIBRARY_SUFFIXES) and ".so." not in base:
                continue
            if base.startswith("minimal_surface."):
                continue  # the extension itself
            # delvewheel: name-<32 hex>.dll; auditwheel: name-<8 hex>.so.N
            base = re.sub(r"-[0-9a-f]{32}(?=\.dll$)", "", base)
            base = re.sub(r"-[0-9a-f]{8}(?=\.so)", "", base)
            names.append(base)
    return sorted(set(names))


def section(title, license_name, source, text):
    rule = "=" * 78
    lines = [rule, title]
    if license_name:
        lines.append("License: " + license_name)
    if source:
        lines.append("Source: " + source)
    lines += [rule, text.strip("\n"), ""]
    return "\n".join(lines) + "\n"


def read_text(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


def eigen_section(eigen_dir):
    """Eigen is compiled into the extension from the copy vendored in the source tree."""
    parts = []
    for name in ("COPYING.MPL2", "COPYING.BSD", "COPYING.README"):
        path = os.path.join(eigen_dir, name)
        if os.path.exists(path):
            parts.append("--- %s ---\n%s" % (name, read_text(path).strip("\n")))
    if not parts:
        raise SystemExit("No Eigen license files found in %s" % eigen_dir)
    parts.append(
        "The Eigen headers used to build this wheel are part of the source distribution of\n"
        "minimal-surface (src/minimal-surface/code/eigen-3.4.0). Eigen's own sources are at\n"
        "https://gitlab.com/libeigen/eigen"
    )
    return section("Eigen 3.4.0 (compiled into the extension)", "MPL-2.0 (some files BSD)",
                   "https://eigen.tuxfamily.org", "\n\n".join(parts))


def conda_packages(prefix):
    metas = []
    meta_dir = os.path.join(prefix, "conda-meta")
    for name in sorted(os.listdir(meta_dir)):
        if name.endswith(".json"):
            with open(os.path.join(meta_dir, name), encoding="utf-8") as f:
                metas.append(json.load(f))
    return metas


def conda_license_section(meta):
    pkg_dir = meta.get("extracted_package_dir")
    lic_dir = os.path.join(pkg_dir, "info", "licenses") if pkg_dir else ""
    if not (lic_dir and os.path.isdir(lic_dir)):
        raise SystemExit("No license files for conda package %s (looked in %s)" % (meta["name"], lic_dir))
    home = ""
    about = os.path.join(pkg_dir, "info", "about.json")
    if os.path.exists(about):
        with open(about, encoding="utf-8") as f:
            home = json.load(f).get("home", "")
    texts = []
    for root, _dirs, files in os.walk(lic_dir):
        for name in sorted(files):
            path = os.path.join(root, name)
            texts.append("--- %s ---\n%s" % (os.path.relpath(path, lic_dir).replace(os.sep, "/"),
                                              read_text(path).strip("\n")))
    if not texts:
        raise SystemExit("Empty license directory for conda package %s" % meta["name"])
    return section("%s %s" % (meta["name"], meta["version"]), meta.get("license", ""), home,
                   "\n\n".join(texts))


def collect_conda(prefix, wheel, eigen_dir):
    metas = conda_packages(prefix)
    owners = {}
    for meta in metas:
        for f in meta.get("files", []):
            base = os.path.basename(f)
            if base.lower().endswith(LIBRARY_SUFFIXES) or ".so." in base:
                owners.setdefault(base.lower(), meta)

    needed = {}
    microsoft = []
    for lib in bundled_libraries(wheel):
        meta = owners.get(lib.lower())
        if meta is not None:
            needed[meta["name"]] = meta
        elif MICROSOFT_RUNTIME.match(lib):
            microsoft.append(lib)
        else:
            raise SystemExit("Cannot tell which conda package provides the bundled library %s" % lib)

    # statically linked into the extension
    for meta in metas:
        if meta["name"] == "libsimpleitk":
            needed[meta["name"]] = meta
    if "libsimpleitk" not in needed:
        raise SystemExit("libsimpleitk is not installed in %s" % prefix)

    parts = []
    for name in sorted(needed):
        parts.append(conda_license_section(needed[name]))
    if microsoft:
        parts.append(section("Microsoft Visual C++ runtime", "", "", MICROSOFT_NOTE + "\nBundled: " + ", ".join(microsoft)))
    parts.append(eigen_section(eigen_dir))
    print("licenses for: " + ", ".join(sorted(needed)) + (", Microsoft runtime" if microsoft else "") + ", Eigen")
    return "".join(parts)


def first_files(root):
    for dirpath, _dirs, files in os.walk(root):
        for name in sorted(files):
            yield os.path.join(dirpath, name)


def collect_linux(wheel, sdk_licenses, eigen_dir):
    parts = []
    for path in sorted(first_files(sdk_licenses)):
        rel = os.path.relpath(path, sdk_licenses).replace(os.sep, "/")
        parts.append(section(rel.split("/")[0] + " (statically linked)", "", "", "--- %s ---\n%s" % (rel, read_text(path))))
    if not parts:
        raise SystemExit("No SDK license files in %s" % sdk_licenses)

    # libraries the repair step copied in come from the system: find the packages that own them
    for lib in bundled_libraries(wheel):
        out = subprocess.run(["ldconfig", "-p"], capture_output=True, text=True, check=True).stdout
        paths = [line.split("=>")[1].strip() for line in out.splitlines() if line.strip().startswith(lib + " ")]
        if not paths:
            raise SystemExit("Cannot find the bundled library %s on this system" % lib)
        real = os.path.realpath(paths[0])
        pkg = subprocess.run(["rpm", "-qf", real], capture_output=True, text=True, check=True).stdout.split()[0]
        files = subprocess.run(["rpm", "-qL", pkg], capture_output=True, text=True, check=True).stdout.split()
        texts = ["--- %s ---\n%s" % (os.path.basename(f), read_text(f).strip("\n")) for f in files if os.path.isfile(f)]
        if not texts:
            raise SystemExit("The package %s that provides %s has no license files" % (pkg, lib))
        parts.append(section("%s (package %s, bundled in the wheel)" % (lib, pkg), "", "", "\n\n".join(texts)))
    parts.append(eigen_section(eigen_dir))
    return "".join(parts)


def record_line(name, data):
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
    return "%s,sha256=%s,%d" % (name, digest, len(data))


def add_to_wheel(wheel, text):
    """Rewrite the wheel with the license file added to its dist-info, in METADATA and in RECORD."""
    tmp = wheel + ".tmp"
    data = text.encode("utf-8")
    with zipfile.ZipFile(wheel) as src:
        infos = src.infolist()
        dist_info = next(i.filename.split("/")[0] for i in infos if i.filename.split("/")[0].endswith(".dist-info"))
        new_name = "%s/licenses/%s" % (dist_info, NAME)
        record_name = dist_info + "/RECORD"
        metadata_name = dist_info + "/METADATA"

        metadata = src.read(metadata_name).decode("utf-8")
        head, sep, body = metadata.partition("\n\n")
        head_lines = head.split("\n")
        if "License-File: " + NAME not in head_lines:
            last = max((i for i, l in enumerate(head_lines) if l.startswith("License-File:")), default=len(head_lines) - 1)
            head_lines.insert(last + 1, "License-File: " + NAME)
        metadata = "\n".join(head_lines) + sep + body
        metadata_data = metadata.encode("utf-8")

        record = [l for l in src.read(record_name).decode("utf-8").splitlines()
                  if l and not l.startswith((metadata_name + ",", new_name + ",", record_name + ","))]
        record += [record_line(metadata_name, metadata_data), record_line(new_name, data), record_name + ",,"]

        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as dst:
            for info in infos:
                if info.filename in (record_name, metadata_name, new_name):
                    continue
                dst.writestr(info, src.read(info.filename), compress_type=zipfile.ZIP_DEFLATED)
            dst.writestr(metadata_name, metadata_data)
            dst.writestr(new_name, data)
            dst.writestr(record_name, ("\n".join(record) + "\n").encode("utf-8"))
    shutil.move(tmp, wheel)
    print("added %s (%d bytes) to %s" % (new_name, len(data), os.path.basename(wheel)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="platform", required=True)
    conda = sub.add_parser("conda")
    conda.add_argument("--prefix", required=True)
    linux = sub.add_parser("linux")
    linux.add_argument("--sdk-licenses", required=True)
    for p in (conda, linux):
        p.add_argument("--wheel", required=True)
        p.add_argument("--eigen", required=True, help="the vendored Eigen directory")
        p.add_argument("--print", action="store_true", help="also print the text")
    args = parser.parse_args(argv)

    if args.platform == "conda":
        text = collect_conda(args.prefix, args.wheel, args.eigen)
    else:
        text = collect_linux(args.wheel, args.sdk_licenses, args.eigen)
    if args.print:
        print(text)
    add_to_wheel(args.wheel, text)


if __name__ == "__main__":
    sys.exit(main())
