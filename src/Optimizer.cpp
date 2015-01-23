#include "Optimizer.h"

Optimizer::Optimizer(ModelDescription& model_, bool drawnPolicy_, double beta_, string folder_): state(model_.getStateVar()), action(model_.getActionVar()), model(model_), folder(folder_), beta(beta_), optAction(model_.getActionVar()), drawnPolicy(drawnPolicy_), nbEdgesModif(0)
{
    gen.seed(time(NULL));
    cout<<"beta = "<<beta<<endl;
}



void Optimizer::computeLandscape(bool all, Value stateObj)
{
    initGraph();


    cout<<"nb edges: "<<num_edges(graph)<<endl;

    Value stateVal(state);
    stateVal.reset();
    if (!all)
    {
        stateVal=stateObj;
        cout<<"Quasi-distance computation for only one objective..."<<endl;
    }

    do
    {

//        if (stateVal[1]==state.cardinality[1]-1)
//            cout<<"\tQuasi-distance: "<<((double)stateVal[0]/state.cardinality(0))*100.0<<"\%"<<endl;

        optimize(stateVal);

        if (!all)
        {
            cout<<"Done."<<endl;
            return;
        }


    } while(stateVal.next());


}


void Optimizer::computeLandscapeT(bool all, Value stateObj)
{
        //Should find a way to construct the single source graph but I don't understand the code...
    initGraph();

    Graph graphT;
    transpose_graph(graph, graphT);

    graph=graphT;
    
    
    cout<<"nb edges: "<<num_edges(graph)<<endl;

    Value stateVal(state);
    stateVal.reset();
    if (!all)
    {
        stateVal=stateObj;
        cout<<"Quasi-distance computation for single source all objectives..."<<endl;
    }

    do
    {

//        if (stateVal[1]==state.cardinality[1]-1)
//            cout<<"\tQuasi-distance: "<<((double)stateVal[0]/state.cardinality(0))*100.0<<"\%"<<endl;

        optimize(stateVal);

        if (!all)
        {
            cout<<"Done."<<endl;
            return;
        }


    } while(stateVal.next());


}

#ifndef _MEMOPTIM
void Optimizer::computePolicy(const Vertex_t& orig, Value& obj)
{
    Value actionVal(action);
    double s;
    std::pair<InEdgeIterator, InEdgeIterator> in_edges_pair;
    Edge_t e;
    map<Value,double>::iterator prob_it;
    Vertex_t dest;
    InEdgeIterator in_i, in_end;
    in_edges_pair=boost::in_edges(orig, graph);
    map<Vertex_t, double>& qdist_obj(getQuasiDist(obj));

    policy[orig]=vector<double>(action.cardinality(),0.0);
    actionVal.reset();
    double sum_p;

    double tmp=0.0;
        //int maxu=0;
    Value maxu(action);
    double maxg=0.0;
    
    
    do
    {
        s=0.0;
        sum_p=0.0;
        for (in_i=in_edges_pair.first; in_i != in_edges_pair.second; ++in_i)
        {

                e = *(in_i);

                if ((prob_it=graph[e].probAction.find(actionVal))!=graph[e].probAction.end())
                {
                    dest = source(e, graph);
                    s+=qdist_obj[dest]* prob_it->second;
                    sum_p+=prob_it->second;

                        //cerr<<prob_it->second<<" ";
                    
                }
                
        }
            //cerr<<endl;
        
            //Test Steve. Semble donner la même chose.
        double g=(model.cost(graph[orig],action)+s-getQuasiDist(obj)[orig]);
        tmp=exp(-beta*g);
        if(tmp>maxg)
        {
            maxg=tmp;
            maxu=actionVal;
        }
        
        
        storeAction(orig, obj, actionVal,s);


    }  while(actionVal.next());

        //FIXME dimension
    bool b=Distributions::normalize(policy[orig]);
        if (!b)
        {
            cout<<"orig : ("<<graph[orig][0]<<" ; "<<graph[orig][1]<<") ; S = "<<s<<" ; Sum P = "<<sum_p<<endl;
        }


            //Could optimize this by doing it in the first loop
            //setMaxPolicy(orig, bestAction(orig));
        setMaxPolicy(orig, maxu);
}


#else
void Optimizer::computePolicy(const Vertex_t& orig, Value& obj)
{}

#endif

void Optimizer::computePolicy(Value& obj)
{

    Value optAction(action);
    Vertex_t orig;

   // vector<double> distrib;




/** WARNING: in_edges instead of out_edges and source instead of target due to graph inversion. **/

    VertexIterator orig_it, last_orig_it;




    for(boost::tie(orig_it, last_orig_it) = boost::vertices(graph); orig_it != last_orig_it; ++orig_it)
    {
        orig=*orig_it;
        if (graph[orig][1]==state.cardinality(1)-1)
            cout<<"\tPolicy: "<<((double)graph[orig][0])/state.cardinality(0)*100.0<<"\%"<<endl;

        computePolicy(orig, obj);

    }


}

void Optimizer::toPlot(Value& obj){
    {
        ofstream file((folder+"/meta.dat").c_str());
            //file<<state.cardinality(0)<<" "<<state.cardinality(1)<<endl;

        file<<state.cardinality(0);
        if(state.nbDim()>1){
            for (int i=1;i<state.nbDim();i++)
                file<<" "<<state.cardinality(i);
        }
        
        file<<endl;
        file<<getType()<<endl;
        file.close();
    }
    {
        ofstream file((folder+"/landscape.dat").c_str());
        map<Vertex_t, double>::iterator it;
        for (it=(getQuasiDist(obj).begin()); it!=(getQuasiDist(obj).end()); it++)
        {
            file<<(*it).second<<" ";
        }
        file.close();
    }
    {
        ofstream file((folder+"/policy.dat").c_str());
        Value s_val(state);
       s_val.reset();
       do
       {
           file<<(getMaxPolicy(getVertex(s_val))[0])<<" ";
       } while(s_val.next());
        file.close();
    }
}


void Optimizer::initGraph()
{
    Value st1(state);
    Value actionVal(action);



    Vertex_t start;

    st1.reset();
    do
    {

        start = boost::add_vertex(graph);
        graph[start].setVariable(state);
        graph[start]=st1;
        setVertex(st1,start);


    } while(st1.next());

    cout<<"nb vertices: "<<num_vertices(graph)<<endl;

    VertexIterator orig, last_orig, dest, last_dest;
    Distributions::Distrib distrib;
    for(boost::tie(orig, last_orig) = boost::vertices(graph); orig != last_orig; ++orig)
    {

        if ((*orig)%(num_vertices(graph)/100)==0)
            cout<<"\tinit: "<<double(*orig)/double(*last_orig)*100.0<<"\%"<<endl;
        actionVal.reset();
        do
        {
            distrib.clear();
            model.getTransDistrib(distrib, graph[*orig],actionVal);
//DEBUG
            
            Distributions::Distrib::iterator tmp_it;
            double sum=0.0;
            
	    
            for(tmp_it=distrib.begin(); tmp_it!=distrib.end(); ++tmp_it){

                    //cerr<<tmp_it->second<<" ";
                sum+=tmp_it->second;
            }
                //cerr<<endl;
            
            ContValue zef(state);
            zef.continuize(graph[*orig]);
            
                
            if(sum<=FLT_EPSILON)
            {
                cerr<<"ERROR "<<graph[*orig][0]<<" "<<graph[*orig][1]<<" "<<graph[*orig][2]<<" "<<zef[0]<<" "<<zef[1]<<" "<<zef[2]<<endl;



            }
                /*
                    //Fix? On reste là où on est? useless    
                Distributions::Distrib d;
                Distributions::dirac(d, state, 0, zef[0]);
                Distributions::dirac(d, state, 1, zef[1]);
                Distributions::dirac(d, state, 2, zef[2]);
                                
                insertDistribInGraph(d, *orig, actionVal);

            }
            
            else
                */
            insertDistribInGraph(distrib, *orig, actionVal);


        } while(actionVal.next()) ;

    }
    cout<<"nb edges: "<<num_edges(graph)<<endl;
//            ofstream file("graph.viz");
//        write_graphviz(file,graph);

}


#ifndef _MEMOPTIM
void Optimizer::removeDistrib(const Vertex_t& orig, Value& actionVal)
{
    std::pair<InEdgeIterator, InEdgeIterator> in_edges_pair;
    InEdgeIterator in_i;
    Edge_t e;
    map<Value,double>::iterator prob_it, it;
    double val;

    in_edges_pair=boost::in_edges(orig, graph);
    for (in_i=in_edges_pair.first; in_i != in_edges_pair.second; ++in_i)
    {
        e = *(in_i);

        if ((prob_it=graph[e].probAction.find(actionVal))!=graph[e].probAction.end())
        {
            val=model.cost(graph[orig],actionVal)/(prob_it->second);
            graph[e].probAction.erase(prob_it);
            if (abs(val-graph[e].value)<FLT_EPSILON)
            {
                if (target(e,graph)==source(e,graph))
                    graph[e].value=0.0;
                else
                    graph[e].value=INF;

                for(it=graph[e].probAction.begin(); it!=graph[e].probAction.end(); ++it)
                {
                    val=model.cost(graph[orig],it->first)/(graph[e].probAction[it->first]);
                    if (val < graph[e].value)
                        graph[e].value=val;

                }
            }

        }
    }
}

#else
void Optimizer::removeDistrib(const Vertex_t& orig, Value& actionVal)
{}
#endif

void Optimizer::insertDistribInGraph(Distributions::Distrib& distrib, const Vertex_t& orig, Value& actionVal)
{

    Edge_t e; bool b;
    OutEdgeIterator out_e_first, out_e_last, out_e_it;
    for(std::map<Value, double>::iterator it=distrib.begin(); it!=distrib.end(); ++it)
    {

        boost::tie(out_e_first,out_e_last)=out_edges(getVertex(it->first),graph);
        bool exist=false;
        for(out_e_it=out_e_first; out_e_it!=out_e_last; ++out_e_it)
        {
            if (target(*out_e_it,graph)==orig)
            {
                exist=true;
                e=*out_e_it;
                break;
            }
        }
        if (!exist)
        {
            boost::tie(e,b) = boost::add_edge(getVertex(it->first),orig,graph);
        }
        initEdge(e, orig, *it, actionVal,exist);

        #ifndef _MEMOPTIM
            //if(it->second<0.001)
            //cerr<<"LOW IT"<<endl;
            //cerr<<it->second<<" ";
        #ifdef _PROBACTIONTHRES
        if(it->second > PROBACTION_THRES)
            #endif
        graph[e].probAction[actionVal]=it->second;
        #endif
    }
}

void OptimizerQuasiDist::initEdge(Edge_t& e, const Vertex_t& orig, std::pair<const Value, double>& p, Value& action, bool& exist)
{

    if (!exist)
    {
        if (target(e,graph)==source(e,graph))
            graph[e].value=0.0;
        else
            graph[e].value=INF;
    }
//                if ( ! p.second>FLT_EPSILON)
//                    cout<<"p petit"<<endl;

    val=model.cost(graph[orig],action)/(p.second);
    //val=model.cost(graph[orig],action)*(1.0-p.second)/(p.second);
    if (val < graph[e].value)
        graph[e].value=val;
    nbEdgesModif++;
}




Optimizer::~Optimizer()
{
    //dtor
}




void OptimizerPureProbs::initEdge(Edge_t& e, VertexIterator& orig, std::pair<const Value, double>& p, Value& action, bool& exist)
{
    val=model.P_UkX(graph[*orig],action)*(p.second);
    if(!exist || val > graph[e].value)
        graph[e].value=val;
}
