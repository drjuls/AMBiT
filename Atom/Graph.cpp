#include "Graph.h"
#include "AvisLib.h"
#include "AvCRT.h"

static IAvGroup* root;
static IAvGraph2D* graph2d;
static bool initialised = false;
static int viewer;
static std::vector<std::wstring> current_datasets;
static std::vector<std::wstring> current_plots;

void Graph::Initialise()
{
    root = (IAvGroup*)avGetObject("/");
    IAvGraphs *graphs;

    root->get_Graphs(&graphs);
    graphs->CreateGraph2D(L"Graph 2D", &graph2d);

    IAvTransform *transform;
    graph2d->get_Transform(&transform);
    transform->put_Origin(AV_XCOORD, -50.);
    transform->put_Origin(AV_YCOORD, -5.);
    transform->put_Size(AV_XCOORD, 1000.);
    transform->put_Size(AV_YCOORD, 10.);

    initialised = true;
}

void Graph::SetData(const std::vector<double>& data, const std::string& name)
{
    if(!initialised)
        Initialise();

    int dim[1];
    dim[0] = data.size();

    double* pData = (double*)avAlloc(1, dim, AV_DOUBLE, name.c_str());
    for(unsigned int i=0; i<dim[0]; i++)
    {   pData[i] = data[i];
    }

    std::wstring wname;
    for(unsigned int i=0; i<name.size(); i++)
        wname.push_back(name[i]);

    std::wstring wstringname(L"/");
    wstringname += wname;

    IAvPlots *plots;
    IAvXYPlot *newplot;
    graph2d->get_Plots(&plots);
    plots->CreatePlot((BSTR)wname.c_str(), AV_XYPLOT, (LPDISPATCH*)&newplot);
    newplot->put_YSource((BSTR)wstringname.c_str());

    IAvTransform *transform;
    IAvPosition *position;
    IAvTransform *graph_transform;

    newplot->get_Transform(&transform);
    graph2d->get_Transform(&graph_transform);
    newplot->get_Position(&position);
    transform->put_Inherit(true);

    double value;
    graph_transform->get_Origin(AV_XCOORD, &value);
    transform->put_Origin(AV_XCOORD, value);
    position->put_Origin(AV_XCOORD, value);
    graph_transform->get_Origin(AV_YCOORD, &value);
    transform->put_Origin(AV_YCOORD, value);
    position->put_Origin(AV_YCOORD, value);

    graph_transform->get_Size(AV_XCOORD, &value);
    transform->put_Size(AV_XCOORD, value);
    position->put_Size(AV_XCOORD, value);
    graph_transform->get_Size(AV_YCOORD, &value);
    transform->put_Size(AV_YCOORD, value);
    position->put_Size(AV_YCOORD, value);
}

void Graph::View()
{
    if(!initialised)
        Initialise();

    viewer = avNewViewer();
    avVisible(viewer, true);
}

void Graph::CloseView()
{
    if(!initialised)
        Initialise();

    avCloseViewer(viewer);
}