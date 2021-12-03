#ifndef SORTED_SET_H
#define SORTED_SET_H


struct ordered_set_node {
	int val;
	int b;                              // balance
	int b_buf;                          // buffer for the balance updates and the direction taken
	struct ordered_set_node * p;        // parent
	struct ordered_set_node * l;        // left
	struct ordered_set_node * r;        // right
};


struct ordered_set {
	long max_size;
	long size;
	struct ordered_set_node * nodes;
	struct ordered_set_node * root;
};


struct ordered_set * ordered_set_new(long max_size);
void ordered_set_destroy(struct ordered_set * t);
void ordered_set_sort(struct ordered_set * t, int * A);
int ordered_set_insert(struct ordered_set * t, int val);


#endif /* SORTED_SET_H */

