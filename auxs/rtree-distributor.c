#include "rtree.h"

void rtree_distributor_rr(rtree_root *root, rtree_node *node) {
	node->server = root->rr_next_server;
	root->rr_next_server++;
	root->rr_next_server = root->rr_next_server % root->servers;
}

void rtree_distributor_space(rtree_root *root, const rtree_node *node) {
}

