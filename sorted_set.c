#include <stdio.h>
#include <stdlib.h>

// #include "debug.h"

#include "sorted_set.h"


struct sorted_set *
sorted_set_new(long max_size)
{
	struct sorted_set * t;
	t = malloc(sizeof(*t));
	t->max_size = max_size;
	t->size = 0;
	t->nodes = malloc(max_size * sizeof(*t->nodes));
	t->root = NULL;
	return t;
}


void
sorted_set_destroy(struct sorted_set * t)
{
	free(t->nodes);
	free(t);
}


// A lot faster than the above, which writes b_buf for every iteration.
static
int
sorted_set_sort_rec(struct sorted_set_node * n, int * A, long i)
{
	if (n->l != NULL)
		i = sorted_set_sort_rec(n->l, A, i);
	A[i++] = n->val;
	if (n->r != NULL)
		i = sorted_set_sort_rec(n->r, A, i);
	return i;
}


void
sorted_set_sort(struct sorted_set * t, int * A)
{
	sorted_set_sort_rec(t->root, A, 0);
}


// n : old root  ,  n->r : new root
static
struct sorted_set_node *
rot_l(struct sorted_set_node * n)
{
	struct sorted_set_node * r, * r_l;
	int b, b_r;
	r = n->r;
	b = n->b;
	b_r = r->b;
	if (b_r < 0)
	{
		n->b = b - 1;                                 //         N(b)                  |                R(nbr)
		if (n->b < 0)                                 //  L               R(b_r)       |       N(nb)               Rr
			r->b = b + b_r - 2;                   //              Rl         Rr    |  L            Rl
		else
			r->b = b_r - 1;
	}
	else
	{
		n->b = b - b_r - 1;
		if (n->b < 0)
			r->b = b - 2;
		else
			r->b = b_r - 1;
	}
	r_l = r->l;
	n->r = r_l;
	if (r_l != NULL)
		r_l->p = n;
	r->l = n;
	n->p = r;
	return r;
}


// n : old root  ,  n->l : new root
static
struct sorted_set_node *
rot_r(struct sorted_set_node * n)
{
	struct sorted_set_node * l, * l_r;
	int b, b_l;
	l = n->l;
	b = n->b;
	b_l = l->b;
	if (b_l < 0)
	{
		n->b = b - b_l + 1;                           //                N(b)            |          L(nbl)
		if (n->b < 0)                                 //       L(b_l)             R     |  Ll                N(nb)
			l->b = b_l + 1;                       //   Ll         Lr                |                Lr         R
		else
			l->b = b + 2;
	}
	else
	{
		n->b = b + 1;
		if (n->b < 0)
			l->b = b_l + 1;
		else
			l->b = b + b_l + 2;
	}
	l_r = l->r;
	n->l = l_r;
	if (l_r != NULL)
		l_r->p = n;
	l->r = n;
	n->p = l;
	return l;
}


static
struct sorted_set_node *
rotate(struct sorted_set_node * n)
{
	struct sorted_set_node * c;
	if (n->b == 2)
	{
		if (n->r->b == -1)
		{
			c = rot_r(n->r);
			n->r = c;
			c->p = n;
			return rot_l(n);
		}
		else
			return rot_l(n);
	}
	else
	{
		if (n->l->b == 1)
		{
			c = rot_l(n->l);
			n->l = c;
			c->p = n;
			return rot_r(n);
		}
		else
			return rot_r(n);
	}
}


int
sorted_set_insert(struct sorted_set * t, int val)
{
	struct sorted_set_node * n, * p, * new;
	long j;

	j = t->size;
	if (j >= t->max_size)
	{
		printf("tree full\n");
		exit(1);
		// error("tree full");
	}

	new = &t->nodes[j];
	new->val = val;
	new->b = 0;
	new->p = NULL;
	new->l = NULL;
	new->r = NULL;

	// Empty tree.
	if (j == 0)
	{
		t->root = new;
		t->size = 1;
		return 1;
	}

	// Insertion can safely fail here (no changes to the structure yet).
	n = t->root;
	p = n;
	while (n != NULL)
	{
		if (val == n->val)   // Value already exists, insert fails.
			return 0;
		p = n;
		if (new->val > n->val)
		{
			n->b_buf = 1;
			n = n->r;
		}
		else
		{
			n->b_buf = -1;
			n = n->l;
		}
	}

	// Insertion was successful.
	n = p;
	new->p = n;
	if (n->b_buf == 1)
		n->r = new;
	else
		n->l = new;
	t->size++;

	p = NULL;
	while (1)
	{
		if (n == NULL)
			return 1;
		n->b += n->b_buf;
		if (n->b == 0)             // tree balance unchanged
			return 1;
		else if (n->b == 2 || n->b == -2)
		{
			p = n->p;
			n = rotate(n);
			break;
		}
		n = n->p;
	}

	// Fix child and parent pointers after rotate().
	n->p = p;
	if (p == NULL)
		t->root = n;
	else if (p->b_buf == 1)
		p->r = n;
	else
		p->l = n;
	return 1;
}

