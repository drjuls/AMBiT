#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>

class Graph
{
public:
    Graph() {}
    ~Graph() {}

    static void SetData(const std::vector<double>& data, const std::string& name);
    static void View();
    static void CloseView();

private:
    static void Initialise();
};

#endif