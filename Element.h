#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>

using namespace std;

// 定义元素及其包含关系的类
class Element {
public:
    string name;
    unordered_set<Element*> containedBy;

    Element() { ; }
    Element(string n) : name(n) {}

    // 递归函数，用于深度优先搜索并将结果存入排序数组
    void dfs(Element* elem, vector<Element*>& sorted, unordered_set<Element*>& visited) {
        if (visited.find(elem) != visited.end()) return;

        visited.insert(elem);

        for (Element* parent : elem->containedBy) {
            dfs(parent, sorted, visited);
        }

        sorted.push_back(elem);
    }

    // 根据包含关系进行拓扑排序
    vector<Element*> topologicalSort(vector<Element*>& elements) {
        vector<Element*> sorted;
        unordered_set<Element*> visited;

        for (Element* elem : elements) {
            dfs(elem, sorted, visited);
        }

        reverse(sorted.begin(), sorted.end());

        return sorted;
    }

};



