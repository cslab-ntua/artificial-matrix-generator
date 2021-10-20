#ifndef SORTED_SET_H
#define SORTED_SET_H


struct sorted_set_node {
	int val;
	int b;                               // balance
	int b_buf;                           // buffer for the balance updates and the direction taken
	struct sorted_set_node * p;        // parent
	struct sorted_set_node * l;        // left
	struct sorted_set_node * r;        // right
};


struct sorted_set {
	long max_size;
	long size;
	struct sorted_set_node * nodes;
	struct sorted_set_node * root;
};


struct sorted_set * sorted_set_new(long max_size);
void sorted_set_destroy(struct sorted_set * t);
void sorted_set_sort(struct sorted_set * t, int * A);
int sorted_set_insert(struct sorted_set * t, int val);


#endif /* SORTED_SET_H */

