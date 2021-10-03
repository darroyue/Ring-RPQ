
/** 
    A class for arrays that can be initialized in constant time (total!)
    Copyright (C) 2021 Diego Arroyuelo

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
**/


#ifndef INIT_ARRAY
#define INIT_ARRAY

#include <bits/stdc++.h> 
#include <cstdint>


using namespace std;

template<typename array_t>
class initializable_array
{
    uint64_t D; // array size
    array_t init_value;  // value used to initialize the array
    
    array_t* V;  // array
    
    uint64_t* S;
    uint64_t top;    
    
    uint64_t* U; // auxiliary array

public:
    
    initializable_array() = default;

    
    initializable_array(const uint64_t _D, array_t _init_value)
    {
        D = _D;
        top = 0;
        V = new array_t[D];
        S = new uint64_t[D];
        U = new uint64_t[D];
        init_value = _init_value;
    };
   
 
    ~initializable_array()
    {
        delete [] V;
        delete [] S;
        delete [] U;
    };


    array_t operator[](const uint64_t i) const
    {
        if (U[i] < top && S[U[i]] == i)
            return V[i];
        else
            return init_value;   
    };


    array_t& operator[](const uint64_t i)   
    {
        if (U[i] < top && S[U[i]] == i)
            return V[i];
        else {
            U[i] = top;
            S[top++] = i;
            return V[i];
        }
    };


    uint64_t size()
    {
        return D;
    };
};
#endif
