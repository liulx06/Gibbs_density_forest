int main(int argc, char** argv);

int BSP_C(vector<string> params);

int BSP_NC(vector<string> params);

int density(vector<string> & params);

int density_2(vector<string> & params);

int density_tree(vector<string> & params);

int countfrompartition(vector<string> params);

int BSP_C_forden(vector<string> params);

void print_usage_and_exit();

int BSP_tree_C(vector<string> params );

int BSP_tree_NC(vector<string> params);

int CDF(vector<string> & params);

int CDF_tree(vector<string> & params);

void test_partition_generation( Partition &mytestpartition ,  Region_Node *  mytestregion_a);

void test_coding( tree	BspTree, string ofilename);

int dec2bin(int n);

