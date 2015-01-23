#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <boost/config.hpp>


//#include <boost/tuple/tuple.hpp>
//#include <boost/tuple/tuple_comparison.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/transpose_graph.hpp>


#include "Specification.h"

#define INF 10e10

#define PROBACTION_THRES 0.001 //threshold for a probaction

using namespace std;
using namespace boost;




//class Vertex
//{
//  public:
//    tuple<int ,int> pos;
//};
//


class Edge
{
    public:

        double value;
        #ifndef _MEMOPTIM
        std::map<Value, double> probAction;  //FIXME optimize memory!
        #endif

        #ifndef _NOSERIAL
    private:
        friend class boost::serialization::access;
        template <class Archive> void serialize(Archive &ar, const unsigned int version)
        {
            #ifndef _MEMOPTIM
            ar & value & probAction; // & policy;
            #else
            ar & value;
            #endif
        }
        #endif
};




class Optimizer
{
    public:
        Optimizer(ModelDescription& model_, bool drawnPolicy_, double beta_=1.0, string folder_=".");
        virtual void computeLandscape(bool all=true, Value stateObj=Value());
        virtual void computeLandscapeT(bool all=true, Value stateObj=Value());
        virtual void computePolicy(Value& obj);

        Value& chooseAction(Value& stateVal, Value& obj)
        {
            if (false) //policyToUpdate[getVertex(stateVal)])
            {
                computePolicy(getVertex(stateVal), obj);
                policyToUpdate[getVertex(stateVal)]=false;
                //cout<<"Update ("<<stateVal[0]<<" ; "<<stateVal[1]<<")"<<endl;
            }

            if (!drawnPolicy)
            {
                    //cout<<"max"<<endl;
                return maxPolicy[getVertex(stateVal)];
            }
            else
            {
                optAction[0]=Distributions::draw(policy[getVertex(stateVal)], gen);
                return optAction;
            }

        }

        void toUpdate()
        {
            for (map<Vertex_t,vector<double> >::iterator it=policy.begin(); it!=policy.end(); ++it)
                policyToUpdate[it->first]=true;
        }

        void load(string fileName)
        {
            cout<<"Loading "<<fileName<<" ..."<<endl;
            ifstream file(fileName.c_str());
            boost::archive::text_iarchive ia(file);

            ia >> *this;
            file.close();
        }

        void save(string fileName)
        {
            ofstream file(fileName.c_str());
            boost::archive::text_oarchive oa(file);

            oa << *this;
            file.close();

            cout<<"Serialized in "<<fileName<<"."<<endl;
        }
        void toPlot(Value& obj);
        void setDir(string dir) {folder=dir;}
        Variable& getStateVar() {return state;}
        Variable& getActionVar() {return action;}
        virtual ~Optimizer();

        typedef boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS, Value, Edge > Graph;
        typedef boost::graph_traits<Graph>::vertex_descriptor Vertex_t;
        typedef boost::graph_traits<Graph>::edge_descriptor Edge_t;
        typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
        typedef boost::graph_traits<Graph>::in_edge_iterator InEdgeIterator;   // because graph is inverted
        typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;

        void computePolicy(const Vertex_t& orig, Value& obj);

        double getBeta() {return beta;}
        void incrBeta(double incr) {beta+=incr;}

            //protected:
        Variable state;
        Variable action;
        ModelDescription& model;
        Graph graph;

      //  vector<map<vector<unsigned int>, vector<double> > > transProbs;
        map<Vertex_t, map<Vertex_t, double> > quasiDist;
        map<vector<unsigned int>, Vertex_t> vertices;
        map<Vertex_t, vector<double> > policy;
        map<Vertex_t, bool > policyToUpdate;
        map<Vertex_t, Value> maxPolicy;
        double val, bestActionScore;

        string folder;

        double beta;

        void initGraph();

        void removeDistrib(const Vertex_t& orig, Value& actionVal);
        void insertDistribInGraph(Distributions::Distrib& distrib, const Vertex_t& orig, Value& actionVal);



        Vertex_t& getVertex(const Value& stateVal)
        {
            return vertices[stateVal.getValueVect()];
        }
        void setVertex(const Value& stateVal, Vertex_t& v)
        {
            vertices[stateVal.getValueVect()]=v;
        }
        map<Vertex_t, double>& getQuasiDist(Value& obj)
        {
            return quasiDist[getVertex(obj)];
        }

        vector<double>& getPolicy(const Vertex_t& state)
        {
            return policy[state];
        }

        Value& getMaxPolicy(Vertex_t& state_)
        {
            return maxPolicy[state_];
        }

        
        void setMaxPolicy(const Vertex_t& state_, Value& action_)
        {
            maxPolicy[state_]=action_;
        }

    private:
        friend class ModelLearner;
        friend class ModelLearnerOneGaussOneDirac;
        virtual void initEdge(Edge_t& e, const Vertex_t& orig, std::pair<const Value, double>& p, Value& action, bool& b)=0;
        virtual void optimize(Value& stateObj)=0;
        
        virtual void storeAction(const Vertex_t& orig, Value& obj, Value& action, double s)=0;
        Value& bestAction(const Vertex_t& orig)
        {
	  bestActionScore=0.0;//-1000000000000000.0; //0.0;
            vector<double>& pol(getPolicy(orig));
	    //why not use max_element?
	    int idx=0;
            for(vector<double>::iterator it=pol.begin(); it!=pol.end(); ++it )
            {
                if (*it > bestActionScore)
                {
		  //optAction[0]=it - pol.begin();//????
		  optAction[0]=idx;// + pol.begin();
                     bestActionScore=*it;
                }
		idx++;
            }
            return optAction;
        }

        virtual unsigned int  getType()=0;  //get type of the subclass (quasidist (0) vs probs (1) )

        friend class boost::serialization::access;
        template <class Archive> void serialize(Archive &ar, const unsigned int version)
        {
            ar & vertices & quasiDist & policy & maxPolicy ; //& graph & quasiDist; // & policy;
        }
        bool b;
        Value optAction;
        boost::mt19937 gen;
        bool drawnPolicy;
    public:
        unsigned int nbEdgesModif;
};


class OptimizerQuasiDist : public Optimizer
{
    public:
        OptimizerQuasiDist(ModelDescription& model, bool drawnPolicy_, double beta_=1.0, string folder_="."): Optimizer(model, drawnPolicy_, beta_, folder_) {;}


    protected:

    private:
        void initEdge(Edge_t& e, const Vertex_t& orig, std::pair<const Value, double>& p, Value& action, bool& b);
        void optimize(Value& stateObj)
        {
                //dijkstra_shortest_paths(graph, getVertex(stateObj), weight_map(get(&Edge::value, graph)).distance_inf(INF).distance_map(associative_property_map<map<Vertex_t,double> >(getQuasiDist(stateObj))));
            dijkstra_shortest_paths(graph, getVertex(stateObj), weight_map(get(&Edge::value, graph)).distance_map(associative_property_map<map<Vertex_t,double> >(getQuasiDist(stateObj))));
            double dist_max=0.0;
            for(map<Vertex_t, double>::iterator it=getQuasiDist(stateObj).begin(); it!=getQuasiDist(stateObj).end(); ++it)
            {
                if(it->second > dist_max)
                    dist_max=it->second;
            }
            distMax=dist_max;
        }
        void storeAction(const Vertex_t& orig, Value& obj, Value& action, double s)
        {
            double g=(model.cost(graph[orig],action)+s-getQuasiDist(obj)[orig]);

                //FIXME: a better warning
                
            //if(g<0.0)
            //    std::cerr<<"G :"<<g<<" "<<s<<std::endl;
                
            //if (g>=0.0)
	      //std::cerr<<"G :"<<g<<" "<<s<<std::endl;

                //Cause problem for infinite dist
                //g=g/distMax;

	    policy[orig][action[0]]= exp(-beta*(g));

	    //else
	      //policy[orig][action[0]]=0.0;
	      //policy[orig][action[0]]=-100000000.0;
        }

        unsigned int  getType() {return 0;}
        double distMax;
};


class DijVis : public default_dijkstra_visitor
{
    public:
    DijVis(Optimizer& E_) : cpt(0), E(E_) {;}
    void examine_edge(Optimizer::Edge_t e, const Optimizer::Graph& G)
    {
  //      if (cpt++>100000)
;//            cout<<"Edge ("<<G[source(e,G)][0]<<","<<G[source(e,G)][1]<<") -> ("<<G[target(e,G)][0]<<","<<G[target(e,G)][1]<<") examined. Target value = "<<E.getQuasiDist(G[target(e,G)])<<endl;
    }
    void edge_relaxed(Optimizer::Edge_t e, const Optimizer::Graph& G)
    {
  //      if (cpt>100000)
 ;//           cout<<"Edge ("<<G[source(e,G)][0]<<","<<G[source(e,G)][1]<<") -> ("<<G[target(e,G)][0]<<","<<G[target(e,G)][1]<<") relaxed. Target value = "<<E.getQuasiDist(G[target(e,G)])<<endl;
    }
    int cpt;
    Optimizer& E;
};

template <class T>
struct times
{
  T operator()(const T& a, const T& b) const {
   return a * b;
  }
};

class OptimizerPureProbs : public Optimizer
{
    public:
        OptimizerPureProbs(ModelDescription& model_, bool drawnPolicy_, double beta_=1.0, string folder_=".") : Optimizer(model_, beta_, drawnPolicy_, folder_) {;}

    protected:

    private:
        void initEdge(Edge_t& e, VertexIterator& orig, std::pair<const Value, double>& p, Value& action, bool& b);
        void optimize(Value& stateObj)
        {
            dijkstra_shortest_paths(graph, getVertex(stateObj), weight_map(get(&Edge::value, graph)).distance_map(associative_property_map<map<Vertex_t,double> >(getQuasiDist(stateObj))).distance_compare(std::greater<double>()).distance_combine(times<double>()).distance_inf(0).distance_zero(1));

        }
        void storeAction(const Vertex_t& orig, Value& obj, Value& action, double s)
        {
            policy[orig][action[0]]=s;
        }
         unsigned int  getType() {return 1;}
};


#endif // Optimizer_H
