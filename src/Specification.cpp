#include "Specification.h"

Variable::Variable(vector<double>& mins_, vector<double>& maxs_, vector<unsigned int>& cardPerDim_, vector<bool>& isCircular_):nbDim_(cardPerDim_.size()), mins(mins_), maxs(maxs_), cardPerDim(cardPerDim_), isCirc(isCircular_)
{
    card=1;
    vector<unsigned int>::iterator it;
    for ( it=cardPerDim.begin() ; it < cardPerDim.end(); it++ )
    {
        card*=*it;
    }

}

Variable Variable::operator^(const Variable &other) const
{
    vector<double> mins_tmp(mins);
    mins_tmp.insert(mins_tmp.end(),other.mins.begin(), other.mins.end());
    vector<double> maxs_tmp(maxs);
    maxs_tmp.insert(maxs_tmp.end(),other.maxs.begin(), other.maxs.end());
    vector<unsigned int> cardPerDim_tmp(cardPerDim);
    cardPerDim_tmp.insert(cardPerDim_tmp.end(),other.cardPerDim.begin(), other.cardPerDim.end());
    vector<bool> isCircular_tmp(isCirc);
    isCircular_tmp.insert(isCircular_tmp.end(),other.isCirc.begin(), other.isCirc.end());

    return Variable(mins_tmp, maxs_tmp, cardPerDim_tmp, isCircular_tmp);

}


bool Value::next()
{
    for (int dim=var->nbDim_-1; dim>=0; dim--)
    {
        if (value[dim]<(var->cardPerDim[dim]-1))
        {
            value[dim]++;
            return true;
        }
        else
        {
            if (dim==0)
                return false;
            else
                value[dim]=0;
        }
    }
    cout<<"Error in Variable::next() !"<<endl;
    return true;
}



void Value::discretize(ContValue& val)
{
    int it;
    double xi;
    bool breaked=false;


    for (unsigned int dim=0; dim<var->nbDim_; dim++)
    {
            
            //if(val.isCircular(dim))
            val.circularize(dim);

        it=0;
        for(xi=var->mins[dim]; xi<(var->maxs[dim]+EPSI); xi+=(var->maxs[dim]-var->mins[dim])/(double)(var->cardPerDim[dim]))
        {
            if( (val[dim]) <= xi){
                if(it==0)
                {
                    value[dim]=it;
                    breaked=true;
                    break;
                }

                else
                {
                    value[dim]=it-1;
                    breaked=true;
                    break;
                }

            }
            it+=1;
        }
            
        if (breaked)
        {
            breaked=false; //raz for other dim
        }
        else
        {
            value[dim]=it-2;
        }
            

    }


}

Value Value::operator^(const Value &other) const
{
    vector<unsigned int> res(value);
    res.insert(res.end(), other.value.begin(), other.value.end());

    Variable* conc= new Variable((*var)^(*(other.var)));
    Value v(*conc);
    v.value=res;
    return v;
}


ContValue ContValue::operator^(const ContValue &other) const
{
    vector<double> res(cvalue);
    res.insert(res.end(), other.cvalue.begin(), other.cvalue.end());

    Variable* conc= new Variable((*var)^(*(other.var)));
    ContValue v(*conc);
    v.cvalue=res;
    return v;
}

void ContValue::continuize(Value& val)
{
    for (unsigned int dim=0; dim<var->nbDim_; dim++)
    {
        cvalue[dim]=var->mins[dim]+(double(val[dim])+0.5)*(var->maxs[dim]-var->mins[dim])/double(var->cardPerDim[dim]);
    }
}

void ContValue::circularize(unsigned int dim)
{
        if (var->isCirc[dim])
        {
            if(cvalue[dim]>var->maxs[dim])
                cvalue[dim]=var->mins[dim]+(cvalue[dim]-var->maxs[dim]);
            else if(cvalue[dim]<var->mins[dim])
                cvalue[dim]=var->maxs[dim]-(var->mins[dim]-cvalue[dim]);
        }
}





void Distributions::gaussian(Distrib1D& P, Variable& X, unsigned int dim, double mu, double sigma, double nbSig_)
{
    if (sigma<FLT_EPSILON)
        return Distributions::dirac(P, X, dim, mu);

    double nbSig;
    //double security_margin=(X.getMaxs(dim)-X.getMins(dim))/(2.0*double(X.cardinality(dim)));  //FLT_EPSILON;

    if ( X.isCircular(dim) && ((2.0*sigma*nbSig_) >= (X.getMaxs(dim)-X.getMins(dim))))
        nbSig=(X.getMaxs(dim)-X.getMins(dim))/(2.0*sigma) - FLT_EPSILON;
    else
        nbSig=nbSig_;


    

    Value X_val(X),tmpd(X);
    ContValue X_cval(X),tmpc(X);
    vector<double> dist(X.cardinality(dim),0.0);
    unsigned int x,left, right;
    double sum=0.0;
    double dist_x_mu, revert_dist_x_mu;
    
    X_cval[dim]=mu-nbSig*sigma;
    X_val.discretize(X_cval);   // TODO: discretize(dim)
    left=X_val[dim];
    X_cval[dim]=mu+nbSig*sigma;
    X_val.discretize(X_cval);   // TODO: discretize(dim)
    right=X_val[dim];
    

            //Fix the buggy behavior of null distrib
    
    if(!(X.isCircular(dim)) && (left==(X.cardinality(dim)-1))) //we are blocked by edges
        return Distributions::dirac(P, X, dim, X.getMins(dim) );
    if(!(X.isCircular(dim)) && (right==0))
        return Distributions::dirac(P, X, dim, X.getMaxs(dim) );


        //WTF? Check that
    if (right+1 != left)
        right=right+1;
    if (right==X.cardinality(dim))
    {
        if (X.isCircular(dim))
            right=0;
        else
             right--;
    }



    for(x=left;  x!=right; X.isCircular(dim) ? x=(x+1)%X.cardinality(dim) : x++)
    {
     //   cout<<left<<" "<<right<<" "<<x<<endl;

        X_val[dim]=x;
        X_cval.continuize(X_val); // TODO: continuize(dim)
        dist_x_mu=abs(X_cval[dim]-mu);
        if (X.isCircular(dim))
        {
            revert_dist_x_mu=abs(X.getMins(dim) - min(X_cval[dim], mu)) + abs(X.getMaxs(dim) - max(X_cval[dim], mu));
            if (revert_dist_x_mu < dist_x_mu)
                dist_x_mu = revert_dist_x_mu;
        }

        dist[x]=(1.0/(sigma*sqrt(2.0*M_PI)))*exp(-((dist_x_mu)*(dist_x_mu))/(2.0*(sigma*sigma)));
        sum+=dist[x];
    }


//    double probMin=1.0/double(10*X.cardinality(dim));
//    double sum_to_remove=0.0;
//    for(x=0; x<X.cardinality(dim); x++)
//    {
//        if ((dist[x]/sum)<=probMin)
//        {
//            sum_to_remove+=dist[x];
//            dist[x]=0.0;
//
//        }
//
//    }
//    sum-=sum_to_remove;

    for(x=left; x!=right; X.isCircular(dim) ? x=(x+1)%X.cardinality(dim) : x++)
    {
      //  if ((dist[x]/sum)>probMin)
            P[x]= dist[x]/sum;
    }
    //cout<<"dim="<<dim<<" ; mu="<<mu<<" ; sig="<<sigma<<endl; //<<" ; X_cval[dim]="<<X_cval[dim]<<" ; X_val[dim]="<<X_val[dim]<<endl;

}

void Distributions::vonMises(Distrib1D& P, Variable& X, unsigned int dim, double mu, double kappa) //, double sigma, double nbSig)
{
    Value X_val(X);
    ContValue X_cval(X);
    vector<double> dist(X.cardinality(dim),0.0);
    unsigned int x;
    double sum=0.0;
    double dist_x_mu;


    for(x=0; x<X.cardinality(dim); x++)
    {

        X_val[dim]=x;
        X_cval.continuize(X_val); // TODO: continuize(dim)
        dist_x_mu=abs(X_cval[dim]-mu);

        dist[x]=exp(kappa*cos(dist_x_mu));

        sum+=dist[x];
    }

    double sum_to_remove=0.0;
    for(x=0; x<X.cardinality(dim); x++)
    {
        if ((dist[x]/sum)<(1.0/double(X.cardinality(dim)*100)))
        {
            sum_to_remove+=dist[x];
            dist[x]=0.0;

        }

    }
    sum-=sum_to_remove;
    for(x=0; x<X.cardinality(dim); x++)
    {
        if ((dist[x]/sum)>=(1.0/double(X.cardinality(dim)*100)))
        {
            P[x]=dist[x]/sum;

        }
    }


}

void Distributions::dirac(Distrib1D& P, Variable& X, unsigned int dim, double val)
{
    Value X_val(X);
    ContValue X_cval(X);
    X_cval[dim]=val;
    X_val.discretize(X_cval);
    P[X_val[dim]]=1.0;

//    cout<<"dim="<<dim<<" ; val="<<val<<" ; X_cval[dim]="<<X_cval[dim]<<" ; X_val[dim]="<<X_val[dim]<<endl;
}

void Distributions::dirac(Distrib& P, Variable& X, unsigned int dim, double val)
{
    Value X_val(X);
    ContValue X_cval(X);
    X_cval[dim]=val;
    X_val.discretize(X_cval);
    P[X_val]=1.0;

//    cout<<"dim="<<dim<<" ; val="<<val<<" ; X_cval[dim]="<<X_cval[dim]<<" ; X_val[dim]="<<X_val[dim]<<endl;
}


void Distributions::uniform(Distrib1D& P, Variable& X, unsigned int dim)
{
    for(unsigned int x=0; x<X.cardinality(dim); x++)
        P[x]=1.0/double(X.cardinality(dim));

//    cout<<"dim="<<dim<<" ; val="<<val<<" ; X_cval[dim]="<<X_cval[dim]<<" ; X_val[dim]="<<X_val[dim]<<endl;
}


Simulator::~Simulator()
{
        //dtor
}


ModelDescription::~ModelDescription()
{
        //dtor
}


