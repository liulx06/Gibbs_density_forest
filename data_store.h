/*
 * File:   data_store.h
 * Author: johnmu
 *
 * Created on June 26, 2012, 2:44 PM
 */

#ifndef DATA_STORE_H
#define	DATA_STORE_H

#include "stl.h"
/*{{{*/
//- by thchiu
/*template <typename T>
class data_store {
private:
    vector<typename vector<T>::const_iterator> data;
public:
    data_store(){

    }

    data_store(vector<T> &d){
        data.reserve(d.size());
        for(typename vector<T>::const_iterator it = d.begin();it != d.end();it++ ){
            data.push_back(it);
        }
    }

    data_store(const data_store<T> &d){
        data.reserve(d.size());

        for(size_t i = 0;i<d.size();i++){
            data.push_back(d.get_ptr(i));
        }
    }

    void operator=(const data_store<T> &d){
        data.clear();
        data.reserve(d.size());

        for(size_t i = 0;i<d.size();i++){
            data.push_back(d.get_ptr(i));
        }
    }

    T operator[](size_t i) const{
        return *(data[i]);
    }

    typename vector<T>::const_iterator get_ptr(size_t i) const{
        return data[i];
    }

    size_t size() const{
        return data.size();
    }


    typename vector<typename vector<T>::const_iterator>::const_iterator begin(){
        return data.begin();
    }

    typename vector<typename vector<T>::const_iterator>::const_iterator end(){
        return data.end();
    }

    void push_back(typename vector<T>::const_iterator a){
        data.push_back(a);
    }

};*//*}}}*/

#define data_store vector   //+ by thchiu


#endif	/* DATA_STORE_H */

