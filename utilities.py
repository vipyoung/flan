"""
Generic utilities.
"""
import time
import datetime as dt
import re
import os
import cPickle as pickle
import geopy
import geopy.distance
from scipy.sparse import find
import networkx as nx
from matplotlib import pyplot as plt
from matplotlib import collections as mc

def geo_distance(pt1, pt2):
	"""
	Compute distance in meters between two gps points.
	:param pt1: LAT, LON (make sure it's not lon, lat)
	:param pt2: Lat, Lon
	:return: distance in meters
	"""
	startpoint = geopy.Point([pt1[1], pt1[0]])
	endpoint = geopy.Point([pt2[1], pt2[0]])
	return geopy.distance.distance(startpoint, endpoint).meters


def save_B(B, fname=None):
	if fname is None:
		return None

	I, J, V = find(B)
	with open(fname, 'w') as g:
		for i in range(len(I)):
			g.write('%s %s %s\n' % (I[i], J[i], V[i]))


def create_result_graph(other, sGs, fname=None):

	rn = nx.DiGraph()
	node2coordinates = {}
	for idx, sg in enumerate(sGs):
		for node in sg.nodes():
			node2coordinates[(idx, node)] = sg.node[node]['coordinates']

	node2index = other['node2index']

	overallnode2localnode = {}
	for (sg_index, localnode), overallnode in node2index.iteritems():
		overallnode2localnode[overallnode] = (sg_index, localnode)

	I, J, V = find(other['B'])

	for i in range(len(I)):
		o_s_n = I[i]
		o_t_n = J[i]
		v = V[i]
		s_n_coor = node2coordinates[overallnode2localnode[o_s_n]]
		t_n_coor = node2coordinates[overallnode2localnode[o_t_n]]
		rn.add_edge(s_n_coor, t_n_coor)

	if fname is not None:
		with open(fname, 'w') as g:
			for s_n_coor, t_n_coor in rn.edges():
				g.write('%s %s %s %s\n' % (s_n_coor[0], s_n_coor[1], t_n_coor[0], t_n_coor[1]))
	return rn

#
# def save_result_graph(other, sGs, fname=None):
# 	if fname is None:
# 		return None
#
# 	node2coordinates = {}
# 	for idx, sg in enumerate(sGs):
# 		for node in sg.nodes():
# 			node2coordinates[(idx, node)] = sg.node[node]['coordinates']
#
# 	node2index = other['node2index']
#
# 	overallnode2localnode = {}
# 	for (sg_index, localnode), overallnode in node2index.iteritems():
# 		overallnode2localnode[overallnode] = (sg_index, localnode)
#
# 	I, J, V = find(other['B'])
# 	with open(fname, 'w') as g:
# 		for i in range(len(I)):
# 			o_s_n = I[i]
# 			o_t_n = J[i]
# 			v = V[i]
# 			s_n_coor = node2coordinates[overallnode2localnode[o_s_n]]
# 			t_n_coor = node2coordinates[overallnode2localnode[o_t_n]]
# 			g.write('%s %s %s %s\n' % (s_n_coor[0], s_n_coor[1], t_n_coor[0], t_n_coor[1]))


def draw_subgraphs(SGs, bbox, showNodes=False):
	fig, ax = plt.subplots(nrows=len(SGs)/3+1, ncols=3, figsize=(30, 40))

	union_lines = []
	for i, G in enumerate(SGs):
		lines = [[s, t] for s, t in G.edges()]
		union_lines += lines
		lc = mc.LineCollection(lines, colors='black', linewidths=1)
		ax1 = plt.subplot(len(SGs)/3+1, 3, i+1)
		ax1.add_collection(lc)
		ax1.autoscale()
		ax1.set_xlim(bbox[1], bbox[0])
		ax1.set_ylim(bbox[3], bbox[2])
		if showNodes == True:
			plt.scatter([node[0] for node in G.nodes()], [node[1] for node in G.nodes()], s=10)

	# Print Union of all subgraphs:
	lc = mc.LineCollection(union_lines, colors='red', linewidths=2)
	ax1 = plt.subplot(len(SGs) / 3 + 1, 3, len(SGs) + 1)
	ax1.add_collection(lc)
	ax1.autoscale()
	ax1.set_xlim(bbox[1], bbox[0])
	ax1.set_ylim(bbox[3], bbox[2])
	if showNodes == True:
		plt.scatter([node[0] for G in SGs for node in G.nodes()], [node[1] for G in SGs for node in G.nodes()], s=10)
	plt.show()


def draw_union_subgraphs(SGs, bbox, showNodes=False):
	fullG = nx.DiGraph()
	for graph_idx, inputG in enumerate(SGs):
		fullG = nx.disjoint_union(fullG, inputG)

	fig, ax = plt.subplots()

	lines = [[fullG.node[s]['coordinates'], fullG.node[t]['coordinates']] for s, t in fullG.edges()]
	lc = mc.LineCollection(lines, colors='black', linewidths=1)
	ax.add_collection(lc)
	ax.autoscale()
	ax.set_xlim(bbox[1], bbox[0])
	ax.set_ylim(bbox[3], bbox[2])
	if showNodes == True:
		plt.scatter([fullG.node[node]['coordinates'][0] for node in fullG.nodes()], [fullG.node[node]['coordinates'][1] for node in fullG.nodes()], s=10)



def draw_graph(G, bbox, showNodes=False):
	fig, ax = plt.subplots()
	lines = [[s, t] for s, t in G.edges()]
	lc = mc.LineCollection(lines, colors='black', linewidths=1)
	ax.add_collection(lc)
	ax.autoscale()
	ax.set_xlim(bbox[1], bbox[0])
	ax.set_ylim(bbox[3], bbox[2])
	if showNodes == True:
		plt.scatter([node[0] for node in G.nodes()], [node[1] for node in G.nodes()], s=10)


def compute_map_similarity_tuples(sGs, save=True, path=None):
	"""
	SOfiane: created this method to compute list of tuple similarities in the format:
										(graph1_index, node1, graph2_index, node2, similarity)
	:param sGs:
	:return:
	"""

	sims = []
	similarity_tuples = []
	for i, g1 in enumerate(sGs):
		for g2 in sGs[i+1:]:
			for node1 in g1.nodes():
				for node2 in g2.nodes():
					dist = geo_distance(node1, node2)
					sim = 1 / float(1 + dist)
					similarity_tuples.append((g1.node[node1]['subgraph'], node1, g2.node[node2]['subgraph'], node2,
					                          dist))
					similarity_tuples.append((g2.node[node2]['subgraph'], node2, g1.node[node1]['subgraph'], node1,
					                          dist))
					sims.append(dist)
	print 'Sofiane: median distance is:', sorted(sims)[len(sims)/2]
	assert (save and not path), "Cannot save the similarity tuples, no path is provided."
	with open(path + 'similarity_tuples.txt', 'w') as fsim:
		for sim_tuple in similarity_tuples:
			fsim.write('%s %s %s %s %s' % sim_tuple)
	return similarity_tuples

def compute_map_similarity_tuples_v1(sGs, path=None, threshold=50):
	"""
	SOfiane: created this method to compute list of tuple similarities in the format:
										(graph1_index, node1, graph2_index, node2, similarity)
	:param sGs:
	:return:
	"""

	sims = []
	similarity_tuples = []
	for i, g1 in enumerate(sGs):
		for g2 in sGs[i+1:]:
			for node1 in g1.nodes():
				# Add identity simiarity
				similarity_tuples.append((g1.node[node1]['subgraph'], node1, g1.node[node1]['subgraph'], node1, 1))
				for node2 in g2.nodes():
					dist = geo_distance(g1.node[node1]['coordinates'], g2.node[node2]['coordinates'])
					if dist > 50:
						continue
					sim = 1 / float(1 + dist)
					similarity_tuples.append((g1.node[node1]['subgraph'], node1, g2.node[node2]['subgraph'], node2,
					                          sim))
					similarity_tuples.append((g2.node[node2]['subgraph'], node2, g1.node[node1]['subgraph'], node1,
					                          sim))
					sims.append(dist)
	print 'Sofiane: median distance is:', sorted(sims)[len(sims)/2]
	assert path is not None, "Cannot save the similarity tuples, no path is provided."
	with open(path + 'similarity_tuples.txt', 'w') as fsim:
		for sim_tuple in similarity_tuples:
			fsim.write('%s %s %s %s %s\n' % sim_tuple)
	return similarity_tuples

def save_sgs_into_edge_file(sgs, path):
	"""
	Sofiane: this is to prepare graphs in the format Eric suggested for his command line code.
	:param sgs:
	:param fname:
	:return:
	"""
	with open(path+'edges.txt', 'w') as fedges,  open(path + 'nodes.txt', 'w') as fnodes:
		for i, sg in enumerate(sgs):
			for s, t in sg.edges():
				fedges.write('%s %s %s 1\n' % (i, s, t))
			for node in sg.nodes():
				fnodes.write('%s %s %s %s\n' % (i, node, sg.node[node]['coordinates'][0], sg.node[node]['coordinates'][1]))



def save_data(data, title='saved_data', date=None):
	date_part = get_date_str(date)
	fname = os.path.join("experiment_results", "{}_{}.pckl".format(title,
																   date_part))
	pickle.dump(data, open(fname, 'wb'))
	print "Wrote the results to: {}".format(fname)
	return fname


def get_date_str(date=None):
	if date is None:
		date = dt.datetime.now()
	date_part = str(date)[:19]
	date_part = re.sub(' ', '_', date_part)
	date_part = re.sub(':', '', date_part)
	return date_part


def cost(assignment, P):
	t0 = time.time()
	X, x = assignment.construct_Xx()
	sum_y = assignment.count_opened_entities()
	unary_sum = (P.D.multiply(X)).sum()
	if P.AA.shape[0] == x.shape[0]:
		binary_sum = P.g * (x.T.dot(P.AA).dot(x))[0, 0] / 2.0
	elif P.A.shape[0] == X.shape[0]:
		qterm = P.A.dot(X).dot(P.A).tolil()
		qterm = qterm.reshape((x.shape[0], 1))
		binary_sum = P.g * (x.T.dot(qterm))[0, 0] / 2.0
	else:
		assert 1==0
	print P.f, sum_y, unary_sum, -binary_sum
	c = P.f*sum_y + unary_sum - binary_sum
	print "Cost took %f seconds." % (time.time()-t0)
	return c


def cost_BROKEN(assignment, P):
	t0 = time.time()
	X, x = assignment.construct_Xx()
	print "Xx", X.shape, x.shape
	sum_y = assignment.count_opened_entities()
	unary_sum = (P.D.multiply(X)).sum()
	assert P.A.shape[0] == X.shape[0]
	qterm = P.A.dot(X).dot(P.A).tolil()
	print x.shape, qterm.shape
	qterm = qterm.reshape((x.shape[0], 1))
	binary_sum = P.g * (x.T.dot(qterm))[0, 0] / 2.0
	print P.f, sum_y, unary_sum, -binary_sum
	c = P.f*sum_y + unary_sum - binary_sum
	print "Cost took %f seconds." % (time.time()-t0)
	return c


def delta_cost(item, dst_clust, assignment, P):
	cur_clust = assignment.matches[item]
	if cur_clust == dst_clust:
		return 0
	c = 0
	if len(assignment.clusters.get(cur_clust, [])) == 1:
		# Emptying cluster
		c -= P.f
	if len(assignment.clusters.get(dst_clust, [])) == 0:
		c += P.f

	# Distance cost
	c += P.D[item,dst_clust] - P.D[item, cur_clust]

	# Binary costs
	for neigh in P.adj_list[item]:
		neigh_clust = assignment.matches[neigh]
		assert P.A[item, neigh] == 1, "Non-neighbor in neighbor list"
		c -= P.g * (P.A[dst_clust, neigh_clust] - P.A[cur_clust, neigh_clust])
	return c


def delta_cost2(item, dst_clust, assignment, P):
	cur_clust = assignment.matches[item]
	if cur_clust == dst_clust:
		return 0
	c = 0
	if len(assignment.clusters.get(cur_clust, [])) == 1:
		# Emptying cluster
		c -= P.f
	if len(assignment.clusters.get(dst_clust, [])) == 0:
		c += P.f

	# Distance cost
	c += P.D[item,dst_clust] - P.D[item, cur_clust]

	# Binary costs
	for other in P.graph2nodes[P.node2graph[item]]:
		if other == item:
			continue
		other_clust = assignment.matches[other]
		if dst_clust == other_clust:
			# Two distinct mapped to the same
			# TODO Should this really cost g and not more (g seems to work
			# better)
			c += P.g
		if other in P.adj_list[item]:
			# Two neighs mapped to non-neighs OR neighs
			c += P.g * (P.A[cur_clust, other_clust] -
						P.A[dst_clust, other_clust])
		elif P.A[dst_clust, other_clust]:
			# Two non-neighs (from the same input graph) mapped to neighs
			c += 0  #P.g# / 5.0
	return c
