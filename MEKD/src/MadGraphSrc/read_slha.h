#ifndef READ_SLHA_H
#define READ_SLHA_H

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <complex>

using namespace std;

class SLHABlock
{
  public:
    SLHABlock(string name = ""){_name = name;}
    ~SLHABlock(){}

    void set_entry(vector<int> indices, complex<double> value);
    complex<double> get_entry(vector<int> indices, complex<double> def_val = 0);
    void set_name(string name) {_name = name;}
    string get_name(){return _name;}
    int get_indices() { return _indices;}

  private:
    string _name;
    map<vector<int>,complex<double> > _entries;
    unsigned int _indices;
};

class SLHAReader
{
  public:
    SLHAReader(string file_name = "")
	{if(file_name != "") read_slha_file(file_name);}

    void read_slha_file(string file_name);
    complex<double> get_block_entry(string block_name, vector<int> indices, 
			   complex<double> def_val = 0);
    complex<double> get_block_entry(string block_name, int index, 
			   complex<double> def_val = 0);
    void set_block_entry(string block_name, vector<int> indices, 
			   complex<double> value);
    void set_block_entry(string block_name, int index, 
			   complex<double> value);
  private:
    map<string, SLHABlock> _blocks;
};


#endif