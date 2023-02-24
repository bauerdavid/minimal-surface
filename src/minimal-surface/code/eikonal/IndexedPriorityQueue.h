#pragma once
//from https://www.geeksforgeeks.org/indexed-priority-queue-with-implementation/
#include <type_traits>
#include <vector>
#include <unordered_map>
#include <limits.h>
#include <stdio.h>
#include <algorithm>
template<typename KeyT, typename ValueT,
    class Hash = std::hash<KeyT>,
    class Comparator = std::less<ValueT>> class IndexedPriorityQueue {
    // Storing indices of values using key
    std::unordered_map<KeyT, size_t, Hash> m;
    public:
        // Container
        std::vector<std::pair<KeyT, ValueT> > v;
    private:
        // Size
        size_t numberOfElement;

        // Creating a instance of Comparator class
        Comparator comp;

        // Max Capacity
        size_t capacity = LLONG_MAX;

        // Obtaing the index value from hash map
        size_t getValueIndex(KeyT key)
        {
            if (m[key] == 0) {
                return -1;
            }
            return v[m[key] - 1];
        }

        // heapify the container
        void heapify(std::vector<std::pair<KeyT, ValueT> >& v,
            long long int heap_size,
            long long index)
        {
            long long leftChild = 2 * index + 1,
                rightChild = 2 * index + 2,
                suitableNode = index;

            if (leftChild < heap_size
                && comp(v[suitableNode].second,
                    v[leftChild].second)) {
                suitableNode = leftChild;
            }

            if (rightChild < heap_size
                && comp(v[suitableNode].second,
                    v[rightChild].second)) {
                suitableNode = rightChild;
            }

            if (suitableNode != index) {

                // swap the value
                std::pair<KeyT, ValueT> temp = v[index];
                v[index] = v[suitableNode];
                v[suitableNode] = temp;

                // updating the map
                m[v[index].first] = index + 1;
                m[v[suitableNode].first]
                    = suitableNode + 1;

                // heapify other affected nodes
                heapify(v, numberOfElement,
                    suitableNode);
            }
        }
    public:

        IndexedPriorityQueue()
        {
            clear();
        }
        template <typename RandomIt>
        IndexedPriorityQueue(RandomIt first, RandomIt last, ValueT default_val) {
            clear();
            std::transform(first, last, back_inserter(v), [default_val](KeyT key) { return std::make_pair(key, default_val); });
            size_t i = 0;
            std::transform(v.begin(), v.end(), inserter(m, m.end()),
                [&i](pair<KeyT, ValueT> kv) {
                    return std::make_pair(kv.first, ++i);
                }
            );
            numberOfElement = v.size();
        }
        template<class KeyIt, class ValueIt>
        IndexedPriorityQueue(KeyIt keys_first, KeyIt keys_last, ValueIt vals_first, bool fix_heap=true) {
            clear();
            std::transform(keys_first, keys_last, vals_first, back_inserter(v), [](KeyT key, ValueT val) { return std::make_pair(key, val); });
            if(fix_heap)
                std::make_heap(v.begin(), v.end(), [this](pair<KeyT, ValueT> l, pair<KeyT, ValueT> r) {return comp(l.second, r.second); });
            size_t i = 0;
            std::transform(v.begin(), v.end(), inserter(m, m.end()),
                [&i](pair<KeyT, ValueT> kv) {
                    return std::make_pair(kv.first, ++i);
                }
            );
            numberOfElement = v.size();
        }
        template<class RandIt>
        IndexedPriorityQueue(RandIt first, RandIt last) {
            clear();
            v.reserve(std::distance(first, last));
            std::copy(first, last, back_inserter(v));
            std::make_heap(v.begin(), v.end(), [this](pair<KeyT, ValueT> l, pair<KeyT, ValueT> r) {return comp(l.second, r.second); });
            size_t i = 0;
            std::transform(v.begin(), v.end(), inserter(m, m.end()),
                [&i](pair<KeyT, ValueT> kv) {
                    return std::make_pair(kv.first, ++i);
                }
            );
            numberOfElement = v.size();
        }

        void clear() {
            numberOfElement = 0;
            m.clear();
            v.clear();
        }

        void push(KeyT key, ValueT value)
        {
            if (numberOfElement == capacity) {
                std::cout << "Overflow";
                return;
            }
            if (contains(key)) {
                std::cout << "Element Already Exists";
                return;
            }

            // Adding element
            v.push_back(std::make_pair(key, value));
            m[key] = ++numberOfElement;

            size_t index = numberOfElement - 1;

            // Comparing to parent node
            while (index != 0
                && comp(v[(index - 1) / 2].second,
                    v[index].second)) {

                // swap the value
                std::pair<KeyT, ValueT> temp = v[index];
                v[index] = v[(index - 1) / 2];
                v[(index - 1) / 2] = temp;

                // updating the map
                m[v[index].first] = index + 1;
                m[v[(index - 1) / 2].first]
                    = (index - 1) / 2 + 1;

                // updating index in map
                index = (index - 1) / 2;
            }
        }

        void pop()
        {
            if (numberOfElement == 0) {
                std::cout << "UnderFlow";
                return;
            }

            // Removing element
            m.erase(v[0].first);
            v[0] = v.back();
            v.pop_back();
            numberOfElement--;
            heapify(v, numberOfElement, 0);
        }

        std::pair<KeyT, ValueT> top() { return v[0]; }

        long long int size() { return numberOfElement; }

        bool empty() { return numberOfElement == 0; }

        void changeAtKey(KeyT key, ValueT value)
        {
            if (!contains(key)) {
                std::cout << "No Such Key Exist";
                return;
            }
            size_t index = m[key] - 1;
            v[index].second = value;

            // Comparing to child nodes
            heapify(v, numberOfElement, index);

            // Comparing to Parent Node
            while (index != 0
                && comp(v[(index - 1) / 2].second,
                    v[index].second)) {

                // swap the value
                std::pair<KeyT, ValueT> temp = v[index];
                v[index] = v[(index - 1) / 2];
                v[(index - 1) / 2] = temp;

                // updating the map
                m[v[index].first] = index + 1;
                m[v[(index - 1) / 2].first]
                    = (index - 1) / 2 + 1;

                // updating index in map
                index = (index - 1) / 2;
            }
        }

        bool contains(KeyT key) {
            return m[key] != 0;
        }

        void push_or_promote(KeyT key, ValueT value) {
            if (contains(key) && comp(v[m[key] - 1].second, value)) {
                changeAtKey(key, value);
            }
            else if (!contains(key))
                push(key, value);
        }

        ValueT& operator[](KeyT key) {
            if (!contains(key))
                throw "Does not contain key";
            return v[m[key] - 1].second;
        }
};