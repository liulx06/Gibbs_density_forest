/*
 *  general_utils.h
 *  seqtools
 *
 *  Created by John Mu on 3/4/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 */
#ifndef GENERAL_UTILS_H

#define GENERAL_UTILS_H

#include "stl.h"

inline void print_data(vector<vector<double> > &data){
    for(int i = 0;i<(int)data.size();i++){
        for(int  j=0;j<(int)data[i].size();j++){
            cerr << data[i][j] << ',';
        }
        cerr << '\n';
    }
}



// Source: http://mlawire.blogspot.com/2009/07/c-whitespace-trimming-functions.html
inline void ltrim(string& str)
{
	string::size_type pos = 0;
	while (pos < str.size() && (isspace(str[pos]))) pos++;
	str.erase(0, pos);
}
inline void rtrim(string& str)
{
	string::size_type pos = str.size();
	while (pos > 0 && (isspace(str[pos - 1]))) pos--;
	str.erase(pos);
}
inline void trim2(string& str)
{
	ltrim(str);
	rtrim(str);
}

inline vector<string> split(const string &s, char delim)
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// Split using all white space as delims
inline vector<string> split(const string &s)
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(ss >> item) {
        elems.push_back(item);
    }
    return elems;
}





#endif




