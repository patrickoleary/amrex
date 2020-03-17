#ifndef AMReX_StateMap_H
#define AMReX_StateMap_H

#include <map>
#include <utility>
#include <string>

namespace amrex {

    // helper to track names and centerings of the avaliable arrays
    class StateMap
    {
    public:
        void Clear(){ this->Map.clear(); }

        int GetIndex(const std::string &arrayName,
           int centering, int &desc, int &comp);

        int Size(int centering){ return this->Map[centering].size(); }

        int GetName(int centering, int id, std::string &name);

        const char *GetAttributesName(int association);

    protected:
        std::map<int,std::map<std::string,std::pair<int,int>>> Map;
    };

}

#endif

