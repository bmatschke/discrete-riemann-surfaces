'''
Basic data structures for discrete Riemann surfaces:
====================================================

A discrete Riemann surface in the sense of this code is a collection of
quadrilaterals (called quads) that are glued together along their edges.
We only consider oriented surfaces, whose graphs are furthermore bipartite.
Furthermore, some algorithms work only when the surface is regular, i.e.
when the intersection of two faces is again a face.

A vertexFlag is a tuple (v,q), where v is one of the four vertices of a quad q.

An edgeFlag is a tuple (e,q), where e is one of the four edges of a quad q.

A vertex of the surface lies in a set of vertexFlags, which we list in ccw order
(starting at the rightmost vertexFlag in case the vertex lies in the interior
of the surface).

An edge of the surface lies in one or two edgeFlags, depending on whether it is
a boundary edge or an interior edge.
'''

import sage.matrix.symplectic_basis;

#precision = 53; #standard precision
precision = 200;
CC = ComplexField(precision);
CIF = ComplexIntervalField(precision);

########################################################################
##### Class Quad #######################################################
########################################################################

class Quad:

	def __init__(self,indexInMesh):
		self.indexInMesh = indexInMesh;

		self.vf0 = VertexFlag(self,0);
		self.vf1 = VertexFlag(self,1);
		self.vf2 = VertexFlag(self,2);
		self.vf3 = VertexFlag(self,3);

		self.ef0 = EdgeFlag(self,0,self.vf1,self.vf0);
		self.ef1 = EdgeFlag(self,1,self.vf1,self.vf2);
		self.ef2 = EdgeFlag(self,2,self.vf3,self.vf2);
		self.ef3 = EdgeFlag(self,3,self.vf3,self.vf0);

		self.vf0.setAdjacentEdgeFlags(self.ef0,self.ef3);
		self.vf1.setAdjacentEdgeFlags(self.ef1,self.ef0);
		self.vf2.setAdjacentEdgeFlags(self.ef2,self.ef1);
		self.vf3.setAdjacentEdgeFlags(self.ef3,self.ef2);

	def setIndexInMesh(self,indexInMesh):
		self.indexInMesh = indexInMesh;

	def index(self):
		return self.indexInMesh;

	def __str__(self):
		return "q"+str(self.index());

	def __repr__(self):
		return self.__str__();

	def vertexFlags(self):
		return [self.vf0,self.vf1,self.vf2,self.vf3];

	def edgeFlags(self):
		return [self.ef0,self.ef1,self.ef2,self.ef3];

	def vertex(self,i):
		return self.vertexFlags()[i].vertex;

	def edge(self,i):
		return self.edgeFlags()[i].edge;

	def commonEdge(self,q):
		for ef in self.edgeFlags():
			efOpp = ef.oppositeEdgeFlag();
			if efOpp != None:
				if efOpp.quad == q:
					return ef.edge;
		return None;

########################################################################
##### Class VertexFlag #################################################
########################################################################

class VertexFlag:
	
	def __init__(self,quad,indexInQuad,efRight=None,efLeft=None):	
		self.quad = quad;
		self.indexInQuad = indexInQuad;
		self.efRight = efRight;
		self.efLeft = efLeft;
		self.vertex = None;

	def setAdjacentEdgeFlags(self,efRight,efLeft):
		self.efRight = efRight;
		self.efLeft = efLeft;

	def index(self):
		return 4*self.quad.index() + self.indexInQuad;

	def color(self):
		return self.indexInQuad % 2;

	def isWhite(self):
		return self.color() == 0;

	def isBlack(self):
		return not self.isWhite();

	def setVertex(self,vertex):
		self.vertex = vertex;

	def __str__(self):
		return "vf"+str(self.index());

	def __repr__(self):
		return self.__str__();

	def theOtherEdgeFlag(self,ef):
		if ef == self.efRight:
			return self.efLeft;
		if ef == self.efLeft:
			return self.efRight;
		raise Exception("The given edge flag "+str(ef)+" is not in the vertex flag "+str(self)+"!");

########################################################################
##### Class EdgeFlag ###################################################
########################################################################

class EdgeFlag:
	
	def __init__(self,quad,indexInQuad,vfBlack,vfWhite):	
		self.quad = quad;
		self.indexInQuad = indexInQuad;
		self.vfBlack = vfBlack;
		self.vfWhite = vfWhite;
		self.edge = None;

	def index(self):
		return 4*self.quad.index() + self.indexInQuad;

	def orientationInQuad(self):
		#We orient the edges from its black to its white endpoint.
		return (-1)^(self.indexInQuad+1);

	def setEdge(self,edge):
		self.edge = edge;
		
	def __str__(self):
		return "ef"+str(self.index());

	def __repr__(self):
		return self.__str__();

	def theOtherVertexFlag(self,vf):
		if self.vfBlack == vf:
			return self.vfWhite;
		elif self.vfWhite == vf:
			return self.vfBlack;
		else:
			raise Exception("The given vertex flag "+str(vf)+" is not in the edge flag "+str(self)+"!");

	def oppositeEdgeFlag(self):
		e = self.edge;
		if e.isBoundaryEdge():
			return None;
		return e.theOtherEdgeFlag(self);

	def vertexFlags(self):
		return [self.vfWhite,self.vfBlack];

	def vfRight(self):
		if self.orientationInQuad() == -1:
			return self.vfBlack;
		else:
			return self.vfWhite;

	def vfLeft(self):
		if self.orientationInQuad() == -1:
			return self.vfWhite;
		else:
			return self.vfBlack;

	def efRight(self):
		return self.vfRight().efRight;

	def efLeft(self):
		return self.vfLeft().efLeft;

	def efOppositeInQuad(self):
		return self.efRight().efRight();

########################################################################
##### Class Vertex #####################################################
########################################################################

class Vertex:

	def __init__(self,vertexFlags=()):
		'''
		We assume that the tuple of vertex flags is given in
		mathematically positive order around the vertex.
		If the vertex is an interior vertex of the mesh,
		then this tuple of vertex flags may start with any of the adjacent flags.
		'''
		
		self.vertexFlags = vertexFlags;
		for vf in self.vertexFlags:
			vf.setVertex(self);

	def setIndexInMesh(self,indexInMesh):
		self.indexInMesh = indexInMesh;

	def index(self):
		return self.indexInMesh;

	def setIndexInMeshAsInteriorVertex(self,indexInMeshAsInteriorVertex):
		self.indexInMeshAsInteriorVertex = indexInMeshAsInteriorVertex;

	def indexAsInteriorVertex(self):
		return self.indexInMeshAsInteriorVertex;

	def setIndexInMeshAsBoundaryVertex(self,indexInMeshAsBoundaryVertex):
		self.indexInMeshAsBoundaryVertex = indexInMeshAsBoundaryVertex;

	def indexAsBoundaryVertex(self):
		return self.indexInMeshAsBoundaryVertex;

	def isInteriorVertex(self):
		'''
		Checks whether the vertex flags close up at the end.
		'''
		
		efR = self.vertexFlags[0].efRight;
		efL = self.vertexFlags[-1].efLeft;
		if efR.edge.isInteriorEdge():
			#Just as a sanity check:
			if efR.edge.theOtherEdgeFlag(efR) != efL:
				raise Exception("Gluing seems erroneous!");
			return True;
		return False;

	def isBoundaryVertex(self):
		return not self.isInteriorVertex();

	def rightEdgeFlagAlongBoundary(self):
		if (not self.isBoundaryVertex()):
			raise Exception(str(self)+" is not a boundary vertex!");
		return self.vertexFlags[0].efRight;

	def leftEdgeFlagAlongBoundary(self):
		if (not self.isBoundaryVertex()):
			raise Exception(str(self)+" is not a boundary vertex!");
		return self.vertexFlags[-1].efLeft;

	def rightAdjacentVertexAlongBoundary(self):
		if (not self.isBoundaryVertex()):
			raise Exception(str(self)+" is not a boundary vertex!");
		vf = self.vertexFlags[0];
		ef = vf.efRight;
		vfNeighbor = ef.theOtherVertexFlag(vf);
		return vfNeighbor;		

	def leftAdjacentVertexAlongBoundary(self):
		if (not self.isBoundaryVertex()):
			raise Exception(str(self)+" is not a boundary vertex!");
		vf = self.vertexFlags[-1];
		ef = vf.efLeft;
		vfNeighbor = ef.theOtherVertexFlag(vf);
		return vfNeighbor;

	def __str__(self):
		return "v"+str(self.index());

	def __repr__(self):
		return self.__str__();

	def adjacentEdges(self):
		'''
		Returns adjacent edges in ccw order.
		'''
		
		result = [vf.efRight.edge for vf in self.vertexFlags];
		if self.isBoundaryVertex():
			result.append(self.leftEdgeFlagAlongBoundary().edge);
		return result;

	def adjacentQuads(self):
		#Returns adjacent quads in ccw order.
		return [vf.quad for vf in self.vertexFlags];

	def isWhite(self):
		return self.vertexFlags[0].isWhite();

	def isBlack(self):
		return self.vertexFlags[0].isBlack();

########################################################################
##### Class Edge #######################################################
########################################################################

class Edge:

	def __init__(self,edgeFlags=()):
		'''
		Make sure that the first edge flag is the one that gives
		a positive boundary orientation to the edge.
		'''
		
		if len(edgeFlags) == 2:
			e1,e2 = edgeFlags;
			if e1.orientationInQuad() == -1:
				e2,e1 = e1,e2;
			if e1.orientationInQuad() == e2.orientationInQuad():
				raise Exception("Glued quads in unoriented fashion!");
			edgeFlags = (e1,e2);			
		self.edgeFlags = edgeFlags;
		for ef in self.edgeFlags:
			ef.setEdge(self);

	def setIndexInMesh(self,indexInMesh):
		self.indexInMesh = indexInMesh;

	def index(self):
		return self.indexInMesh;

	def setIndexInMeshAsInteriorEdge(self,indexInMeshAsInteriorEdge):
		self.indexInMeshAsInteriorEdge = indexInMeshAsInteriorEdge;

	def indexAsInteriorEdge(self):
		return self.indexInMeshAsInteriorEdge;

	def setIndexInMeshAsBoundaryEdge(self,indexInMeshAsBoundaryEdge):
		self.indexInMeshAsBoundaryEdge = indexInMeshAsBoundaryEdge;

	def indexAsBoundaryEdge(self):
		return self.indexInMeshAsBoundaryEdge;

	def leftEdgeFlag(self):
		if (not self.isInteriorEdge()):
			raise Exception("Edge is assumed to lie in two quads!");
		return self.edgeFlags[0];

	def rightEdgeFlag(self):
		if (not self.isInteriorEdge()):
			raise Exception("Edge is assumed to lie in two quads!");
		return self.edgeFlags[1];

	def leftQuad(self):
		return self.leftEdgeFlag().quad;

	def rightQuad(self):
		return self.rightEdgeFlag().quad;

	def isInteriorEdge(self):
		return len(self.edgeFlags) == 2;

	def isBoundaryEdge(self):
		return not self.isInteriorEdge();

	def whiteVertex(self):
		return self.edgeFlags[0].vfWhite.vertex;

	def blackVertex(self):
		return self.edgeFlags[0].vfBlack.vertex;

	def vertices(self):
		return [self.whiteVertex(),self.blackVertex()];

	def __str__(self):
		return "e"+str(self.index());

	def __repr__(self):
		return self.__str__();

	def theOtherVertex(self,v):
		if self.whiteVertex() == v:
			return self.blackVertex();
		elif self.blackVertex() == v:
			return self.whiteVertex();
		else:
			raise Exception("The given vertex "+str(v)+" is not in the edge "+str(self)+"!");

	def theOtherEdgeFlag(self,ef):
		if (not self.isInteriorEdge()):
			raise Exception("Edge is assumed to lie in two quads!");
		if self.edgeFlags[0] == ef:
			return self.edgeFlags[1];
		elif self.edgeFlags[1] == ef:
			return self.edgeFlags[0];
		else:
			raise Exception("The given edge flag "+str(ef)+" is not in the edge "+str(self)+"!");

	def vertexFlagTowardsOtherEdge(self,e):
		for ef in self.edgeFlags:
			if ef.vfWhite.efLeft.edge == e or ef.vfWhite.efRight.edge == e:
				return ef.vfWhite;
			if ef.vfBlack.efLeft.edge == e or ef.vfBlack.efRight.edge == e:
				return ef.vfBlack;
		return None;

	def adjacentEdges(self):
		edgesWhite = self.whiteVertex().adjacentEdges();
		edgesWhite.remove(self);
		edgesBlack = self.blackVertex().adjacentEdges();
		edgesBlack.remove(self);
		return edgesWhite + edgesBlack;

	def vertexFlags(self):
		return [vf for ef in self.edgeFlags for vf in ef.vertexFlags()];

########################################################################
##### Class DiscreteRiemannSurface #####################################
########################################################################

class DiscreteRiemannSurface:

	def __init__(self,numQuads=1,pairsOfIndicesOfGluedEdgeFlags=[]):
		'''
		The quads will be numbered from 0 up to (numQuads-1).
		The edge flags will be numbered from 0 up to 4*(numQuads-1),
		such that the first four edges belong to the first quad, etc.
		'''
		
		self.quads = [];
		self.edgeFlags = [];
		self.vertexFlags = [];
		for i in range(numQuads):
			q = Quad(i);
			self.quads.append(q);
			self.edgeFlags += q.edgeFlags();
			self.vertexFlags += q.vertexFlags();

		self.edges = [Edge(edgeFlags=(ef,)) for ef in self.edgeFlags];
		self.vertices = [Vertex(vertexFlags=(vf,)) for vf in self.vertexFlags];
		
		for ef1i,ef2i in pairsOfIndicesOfGluedEdgeFlags:
			ef1 = self.edgeFlags[ef1i];
			ef2 = self.edgeFlags[ef2i];
			self.glueEdgeFlags(ef1,ef2,updateIndicesAtTheEnd=False);
		self.updateIndices();

	def __copy__(self):
		'''
		Warning: This yields not an exact copy! Indices of edges and vertices might be different!
		'''
		
		numQuads = self.f2();
		pairsOfIndicesOfGluedEdgeFlags = [(e.rightEdgeFlag().index(),e.leftEdgeFlag().index()) for e in self.edges if e.isInteriorEdge()];
		return DiscreteRiemannSurface(numQuads,pairsOfIndicesOfGluedEdgeFlags);

	def f0(self):
		return len(self.vertices);

	def f0i(self):
		return len(self.interiorVertices);

	def f0b(self):
		return len(self.boundaryVertices);

	def f02(self):
		return len(self.vertexFlags);

	def f1(self):
		return len(self.edges);

	def f1i(self):
		return len(self.interiorEdges);

	def f1b(self):
		return len(self.boundaryEdges);

	def f12(self):
		return len(self.edgeFlags);
	
	def f2(self):
		return len(self.quads);

	def vertexDegrees_edges(self):
		return [len(v.adjacentEdges()) for v in self.vertices];

	def vertexDegrees_quads(self):
		return [len(v.adjacentQuads()) for v in self.vertices];

	def eulerCharacteristic(self):
		return self.f0() - self.f1() + self.f2();

	def __str__(self):
		return "Discrete Riemann surface with f = ("+str(self.f0())+","+str(self.f1())+","+str(self.f2())+"), chi = "+str(self.eulerCharacteristic())+", g = "+str(self.genus())+" and b = "+str(self.numBoundaryCycles())+".";

	def __repr__(self):
		return self.__str__();

	def isRegular(self):
		'''
		Checks whether each quad has four distinct vertices and that
		the edges are determined by their vertices.
		'''
		
		#Check that each quad has four different vertices:
		for q in self.quads:
			v0 = q.vf0.vertex;
			v1 = q.vf1.vertex;
			v2 = q.vf2.vertex;
			v3 = q.vf3.vertex;
			if v0 == v2 or v1 == v3:
				return False;
		#Check that edges are determined by their vertices:
		for e1 in self.edges:
			for e2 in self.edges:
				if e1.index() < e2.index():
					if e1.whiteVertex() == e2.whiteVertex():
						if e1.blackVertex() == e2.blackVertex():
							return False;
		return True;		

	def computeBoundary(self):
		self.boundaryVertices = [v for v in self.vertices if v.isBoundaryVertex()];
		self.boundaryEdges = [e for e in self.edges if e.isBoundaryEdge()];
		boundaryGraphEdges = [(e.blackVertex(),e.whiteVertex()) for e in self.boundaryEdges];
		self.boundaryGraph = Graph([self.boundaryVertices,boundaryGraphEdges]);

	def numBoundaryCycles(self):
		return self.boundaryGraph.connected_components_number();

	def genus(self):
		'''
		g = (2 - chi - b)/2, i.e.
		2 - 2*g - b = chi
		'''
		return ZZ((2-self.eulerCharacteristic()-self.numBoundaryCycles())/2);

	def updateIndices(self):
		#TODO: Should first order quads, edges and vertices canonically,
		#such that self.copy() yields a copy with the same indices.
		
		for i in range(len(self.quads)):
			self.quads[i].setIndexInMesh(i);
		for i in range(len(self.edges)):
			self.edges[i].setIndexInMesh(i);
		for i in range(len(self.vertices)):
			self.vertices[i].setIndexInMesh(i);
		#Indices of self.edgeFlags are automatically computed via self.quads.
		#Indices of self.vertexFlags are automatically computed via self.quads.

		self.computeBoundary();
		self.interiorEdges = [e for e in self.edges if e.isInteriorEdge()];
		for i in range(len(self.interiorEdges)):
			self.interiorEdges[i].setIndexInMeshAsInteriorEdge(i);
		self.boundaryEdges = [e for e in self.edges if e.isBoundaryEdge()];
		for i in range(len(self.boundaryEdges)):
			self.boundaryEdges[i].setIndexInMeshAsBoundaryEdge(i);

		self.interiorVertices = [v for v in self.vertices if v.isInteriorVertex()];
		for i in range(len(self.interiorVertices)):
			self.interiorVertices[i].setIndexInMeshAsInteriorVertex(i);
		self.boundaryVertices = [e for e in self.vertices if e.isBoundaryVertex()];
		for i in range(len(self.boundaryVertices)):
			self.boundaryVertices[i].setIndexInMeshAsBoundaryVertex(i);

		self.H1_lambda = None;
		self.H1_lambda_matrix = None;
		self.H1_CX = None;
		self.H1_CX_matrix = None;
		self.H1_X = None;
		self.H1_X_matrix = None;
		self.Coh1_lambda = None;
		self.Coh1_lambda_matrix = None;
		self.Coh1_CX = None;
		self.Coh1_CX_matrix = None;
		self.Coh1_X = None;
		self.Coh1_X_matrix = None;

	def glueEdgeFlags(self,ef1,ef2,updateIndicesAtTheEnd=True):
		if ef1.orientationInQuad() == ef2.orientationInQuad() == 0:
			raise Exception("Need to keep the graph bipartite!");

		if ef1.orientationInQuad() == -1:
			ef1,ef2 = ef2,ef1;

		#print "Edge classes:",self.edgeClasses;

		#Update edge identifications and equivalence classes:
		self.edges.remove(ef1.edge);
		self.edges.remove(ef2.edge);
		self.edges.append(Edge(edgeFlags=(ef1,ef2)));

		#print "Vertex classes:",self.vertexClasses;

		#Update vertex equivalence classes:
		vfb1 = ef1.vfBlack;
		vfw1 = ef1.vfWhite;
		vfb2 = ef2.vfBlack;
		vfw2 = ef2.vfWhite;
		if vfb1.vertex != vfb2.vertex: #I.e. the gluing didn't close up the cone 
			self.vertices.remove(vfb1.vertex);
			self.vertices.remove(vfb2.vertex);
			self.vertices.append(Vertex(vfb2.vertex.vertexFlags+vfb1.vertex.vertexFlags)); #Order is important!		
		if vfw1.vertex != vfw2.vertex: #I.e. the gluing didn't close up the cone 
			self.vertices.remove(vfw1.vertex);
			self.vertices.remove(vfw2.vertex);
			self.vertices.append(Vertex(vfw1.vertex.vertexFlags+vfw2.vertex.vertexFlags)); #Order is important!		
		#Update indices of edges and vertices:
		if updateIndicesAtTheEnd:
			self.updateIndices();

	def attach(self,P2,pairsOfIndicesOfGluedEdgeFlags=[]):
		P2 = copy(P2);
		self.quads += P2.quads;
		self.edgeFlags += P2.edgeFlags;
		self.vertexFlags += P2.vertexFlags;
		self.edges += P2.edges;
		self.vertices += P2.vertices;
		
		for ef1i,ef2i in pairsOfIndicesOfGluedEdgeFlags:
			ef1 = self.edgeFlags[ef1i];
			ef2 = P2.edgeFlags[ef2i];
			self.glueEdgeFlags(ef1,ef2,updateIndicesAtTheEnd=False);
		self.updateIndices();

	def subdivision(self,n=2):
		'''
		Returns a surface which is given by subdividing each quad in self
		into n times n pieces.
		'''
		
		gluing = [];
		f2 = self.f2();
		newf2 = n*n*f2;
		#Glue n*n quads together for each quad in self:
		ei = {0:[0,1,2,3], 1:[3,0,1,2]};
		for qi in range(f2):
			q0 = n*n*qi;
			e0 = 4*q0;
			for x in range(n):
				x2 = x % 2;
				for y in range(n):
					y2 = y % 2;
					xy2 = (x+y) % 2;
					xy2o = (x+y+1) % 2;
					if x<n-1:
						gluing.append((e0+4*(x+n*y)+ei[xy2][1],e0+4*(x+1+n*y)+ei[xy2o][3]));
					if y<n-1:
						gluing.append((e0+4*(x+n*y)+ei[xy2][2],e0+4*(x+n*y+n)+ei[xy2o][0]));
		#Subdivide edges of self:
		indices = {};
		indices[0] = [4*x+ei[x%2][0] for x in range(n)];
		indices[1] = [4*(y*n+(n-1))+ei[(n-1+y)%2][1] for y in range(n)];
		indices[2] = [4*(x+n*(n-1))+ei[(n-1+x)%2][2] for x in range(n-1,-1,-1)];
		indices[3] = [4*y*n+ei[y%2][3] for y in range(n-1,-1,-1)];
		for e in self.interiorEdges:
			ef0, ef1 = e.edgeFlags;
			q0i = ef0.quad.index();
			q1i = ef1.quad.index();
			for x in range(n):
				gluing.append((4*n^2*q0i+indices[ef0.indexInQuad][x],4*n^2*q1i+indices[ef1.indexInQuad][n-1-x]));
		#print gluing;
		result = DiscreteRiemannSurface(numQuads=newf2,pairsOfIndicesOfGluedEdgeFlags=gluing);
		return result;

	def DEPRECIATED_subdivision(self):
		'''
		Returns a surface which is given by subdividing each quad in self
		into four pieces.
		'''
		
		gluing = [];
		f2 = self.f2();
		newf2 = 4*f2;
		#Glue 4 quads together for each quad in self:
		for i in range(f2):
			gluing.append((16*i+ 1,16*i+ 6));
			gluing.append((16*i+ 5,16*i+ 8));
			gluing.append((16*i+11,16*i+12));
			gluing.append((16*i+ 2,16*i+15));
		#Subdivide edges of self:
		leftIndex = [7,9,13,3];
		rightIndex = [0,4,10,14];
		for e in self.interiorEdges:
			ef0, ef1 = e.edgeFlags;
			q0i = ef0.quad.index();
			q1i = ef1.quad.index();
			gluing.append((16*q0i+leftIndex[ef0.indexInQuad],16*q1i+rightIndex[ef1.indexInQuad]));
			gluing.append((16*q0i+rightIndex[ef0.indexInQuad],16*q1i+leftIndex[ef1.indexInQuad]));
		print gluing;
		result = DiscreteRiemannSurface(numQuads=newf2,pairsOfIndicesOfGluedEdgeFlags=gluing);
		return result;
	
	def graph_lambda(self):
		'''
		The 1-skeleton of the surface.
		'''
		
		nodes = copy(self.vertices);
		edges = [(e.blackVertex(),e.whiteVertex()) for e in self.edges]
		#print "Nodes:",nodes;
		#print "Edges:",edges;
		G = Graph([nodes,edges]);
		return G;

	def dualGraph_lambda(self):
		nodes = copy(self.quads);
		edges = [(e.leftQuad(),e.rightQuad()) for e in self.edges if e.isInteriorEdge()];
		Gdual = Graph([nodes,edges]);
		return Gdual;

	def graph_X(self):
		'''
		Edge graph.
		'''
		
		nodes = copy(self.edges);
		edges = [(vf.efRight.edge,vf.efLeft.edge) for vf in self.vertexFlags];
		#print "Nodes:",nodes;
		#print "Edges:",edges;
		G = Graph([nodes,edges]);
		return G;

	def graph_CX(self):
		'''
		edgeFlag graph:
		'''
				
		nodes = copy(self.edgeFlags);
		edges = [(vf.efRight,vf.efLeft) for vf in self.vertexFlags];
		edges += [e.edgeFlags for e in self.edges if e.isInteriorEdge()];
		#print "Nodes:",nodes;
		#print "Edges:",edges;
		G = Graph([nodes,edges]);
		return G;

	def boundaryOperatorMatrix_1_to_0_lambda(self):
		M = zero_matrix(ZZ,self.f0(),self.f1());
		for e in self.edges:
			M[e.blackVertex().index(),e.index()] -= 1;
			M[e.whiteVertex().index(),e.index()] += 1;
		return M;

	def boundaryOperatorMatrix_2_to_1_lambda(self):
		M = zero_matrix(ZZ,self.f1(),self.f2());
		for q in self.quads:
			M[q.edge(0).index(),q.index()] -= 1;
			M[q.edge(1).index(),q.index()] += 1;
			M[q.edge(2).index(),q.index()] -= 1;
			M[q.edge(3).index(),q.index()] += 1;
		return M;

	def chainComplex_lambda(self):
		D10 = self.boundaryOperatorMatrix_1_to_0_lambda();
		D21 = self.boundaryOperatorMatrix_2_to_1_lambda();
		#print "D10:",D10;
		#print "D21:",D21;
		#print (D10*D21).is_zero();
		C = ChainComplex(data={2:D21,1:D10},degree=-1);
		return C;

	def homology_lambda(self,deg=None):
		return self.chainComplex_lambda().homology(generators=True,deg=deg);

	def homologyBasis_lambda(self):
		H1 = self.homology_lambda(deg=1);
		basis = [c[1].vector(1) for c in H1];
		B = matrix(basis);
		B = B.LLL();
		return B.rows();

	def intersectionNumber_lambda(self,gamma1,gamma2):
		'''
		Computes the intersection number (in ZZ) between
		two one-cycles gamma1 and gamma2 in the graph (lambda) of the surface.
		'''
		
		result = 0;
		for v in self.vertices:
			sum1 = 0;
			if v.isWhite():
				for e in v.adjacentEdges():
					sum1 += gamma1[e.index()];
					result += sum1 * gamma2[e.index()];
			else:
				for e in v.adjacentEdges():
					result += sum1 * gamma2[e.index()];
					sum1 += gamma1[e.index()];
		return result;

	def intersectionForm_lambda(self,basis=None):
		if basis == None:
			basis = self.homologyBasis_lambda();
		M = matrix(ZZ,len(basis),len(basis));
		for i in range(len(basis)):
			for j in range(len(basis)):
				M[i,j] = self.intersectionNumber_lambda(basis[i],basis[j]);
		return M;

	def normalizedHomologyBasis_lambda(self,basis1=None):
		if basis1 == None:
			basis1 = self.homologyBasis_lambda();
		B1 = matrix(basis1);
		M = self.intersectionForm_lambda(basis1);
		Mnew,C = sage.matrix.symplectic_basis.symplectic_basis_over_ZZ(M);
		#Here Mnew = C*M*C^t.
		return (C*B1).rows();		

	def cohomologyBasis_lambda(self,homologyBasis=None):
		'''
		If a homologyBasis is given, we compute a cohomology basis
		which is dual to it.
		'''
		
		Coh1 = self.chainComplex_lambda().dual().homology(generators=True,deg=1);
		cohBasis = [a[1].vector(1) for a in Coh1];
		cohB = matrix(cohBasis);
		if homologyBasis == None:
			cohB = cohB.LLL();
		else:
			B = matrix(homologyBasis);
			T = cohB * B.transpose();
			cohB = T.inverse() * cohB;
		return cohB.rows();

	def boundaryOperatorMatrix_1_to_0_CX(self):
		f02 = self.f02();
		f1i = self.f1i();
		M = zero_matrix(ZZ,self.f12(),f02+f1i);
		for vf in self.vertexFlags:
			if vf.isWhite():
				M[vf.efLeft.index(),vf.index()] -= 1;
				M[vf.efRight.index(),vf.index()] += 1;
			else:
				M[vf.efLeft.index(),vf.index()] += 1;
				M[vf.efRight.index(),vf.index()] -= 1;
		for e in self.interiorEdges:
			M[e.leftEdgeFlag().index(),f02+e.indexAsInteriorEdge()] -= 1;
			M[e.rightEdgeFlag().index(),f02+e.indexAsInteriorEdge()] += 1;
		return M;

	def boundaryOperatorMatrix_2_to_1_CX(self):
		f02 = self.f02();
		f1i = self.f1i();
		f2 = self.f2();
		f0i = self.f0i();
		M = zero_matrix(ZZ,f02+f1i,f2+f0i);
		for q in self.quads:
			M[q.vf0.index(),q.index()] += 1;
			M[q.vf1.index(),q.index()] -= 1;
			M[q.vf2.index(),q.index()] += 1;
			M[q.vf3.index(),q.index()] -= 1;
		for v in self.interiorVertices:
			for vf in v.vertexFlags:
				if vf.isWhite():
					M[vf.index(),f2+v.indexAsInteriorVertex()] -= 1;
					M[f02+vf.efLeft.edge.indexAsInteriorEdge(),f2+v.indexAsInteriorVertex()] +=1;
				else:
					M[vf.index(),f2+v.indexAsInteriorVertex()] += 1;
					M[f02+vf.efLeft.edge.indexAsInteriorEdge(),f2+v.indexAsInteriorVertex()] -= 1;
		return M;

	def chainComplex_CX(self):
		D10 = self.boundaryOperatorMatrix_1_to_0_CX();
		D21 = self.boundaryOperatorMatrix_2_to_1_CX();
		C = ChainComplex(data={2:D21,1:D10},degree=-1);
		return C;

	def boundaryOperatorMatrix_1_to_0_X(self):
		f02 = self.f02();
		M = zero_matrix(ZZ,self.f12(),f02);
		for vf in self.vertexFlags:
			if vf.isWhite():
				M[vf.efLeft.index(),vf.index()] -= 1;
				M[vf.efRight.index(),vf.index()] += 1;
			else:
				M[vf.efLeft.index(),vf.index()] += 1;
				M[vf.efRight.index(),vf.index()] -= 1;
		return M;

	def boundaryOperatorMatrix_2_to_1_X(self):
		f02 = self.f02();
		f2 = self.f2();
		f0i = self.f0i();
		M = zero_matrix(ZZ,f02,f2+f0i);
		for q in self.quads:
			M[q.vf0.index(),q.index()] += 1;
			M[q.vf1.index(),q.index()] -= 1;
			M[q.vf2.index(),q.index()] += 1;
			M[q.vf3.index(),q.index()] -= 1;
		for v in self.interiorVertices:
			for vf in v.vertexFlags:
				if vf.isWhite():
					M[vf.index(),f2+v.indexAsInteriorVertex()] -= 1;
				else:
					M[vf.index(),f2+v.indexAsInteriorVertex()] += 1;
		return M;

	def chainComplex_X(self):
		D10 = self.boundaryOperatorMatrix_1_to_0_X();
		D21 = self.boundaryOperatorMatrix_2_to_1_X();
		C = ChainComplex(data={2:D21,1:D10},degree=-1);
		return C;

	def intersectionNumber_X(self,gamma1,gamma2):
		'''
		Computes the intersection number (in ZZ) between
		two one-cycles gamma1 and gamma2 in the edge graph X of the surface.
		'''
		
		result = 0;
		for e in self.edges:
			sum1 = 0;
			for ef in e.edgeFlags:
				sum1 -= ef.orientationInQuad() * gamma1[ef.vfBlack.index()];
			sum2 = 0;
			ef = e.edgeFlags[0];
			sum2 += ef.orientationInQuad() * gamma2[ef.vfWhite.index()];
			sum2 += ef.orientationInQuad() * gamma2[ef.vfBlack.index()];
			result += sum1 * sum2;
		return result;

	def intersectionForm_X(self,basis=None):
		if basis == None:
			basis = self.H1_X;
		M = matrix(ZZ,len(basis),len(basis));
		for i in range(len(basis)):
			for j in range(len(basis)):
				M[i,j] = self.intersectionNumber_X(basis[i],basis[j]);
		return M;

	def computeHomologies(self):
		#Compute H_1(surface) e wrt. the celldecomposition given by lambda:
		C_lambda = self.chainComplex_lambda();
		H1_lambda = [c[1].vector(1) for c in C_lambda.homology(generators=True,deg=1)];
		H1_lambda = matrix(H1_lambda).LLL().rows();
		#Normalize first homology basis:
		I_lambda = self.intersectionForm_lambda(basis=H1_lambda);
		I_normalized_lambda,Trf = sage.matrix.symplectic_basis.symplectic_basis_over_ZZ(I_lambda);
		#Here, I_normalized_lambda = Trf * I_lambda * Trf^{-1} is the representing matrix of the standard symplectic form.
		H1_lambda_matrix = Trf * matrix(H1_lambda);
		H1_lambda = H1_lambda_matrix.rows();
		self.H1_lambda_matrix = H1_lambda_matrix;
		self.H1_lambda = H1_lambda;

		#print 1;
		
		#Compute H^1(surface) wrt. the celldecomposition given by lambda:
		Coh1_lambda = [a[1].vector(1) for a in C_lambda.dual().homology(generators=True,deg=1)];
		Coh1_lambda_matrix = matrix(Coh1_lambda);
		#Choose basis that makes it dual to the chosen homology basis:
		Trf = Coh1_lambda_matrix * H1_lambda_matrix.transpose();
		Coh1_lambda_matrix = Trf.inverse() * Coh1_lambda_matrix;
		Coh1_lambda = Coh1_lambda_matrix.rows();
		self.Coh1_lambda_matrix = Coh1_lambda_matrix;
		self.Coh1_lambda = Coh1_lambda;

		#print 2;

		#Compute H^1(surface) wrt. CX, the edgeFlag-graph:
		#The edges in CX are first indexed by the vertexFlags and then by the edges.
		#We use a canonical map
		f02 = self.f02();
		f1i = self.f1i();
		Coh1_CX = [];
		for a_lambda in Coh1_lambda:
			a_CX = zero_vector(f02+f1i);
			for vf in self.vertexFlags:
				a_CX[vf.index()] = a_lambda[vf.efLeft.edge.index()];
			for e in self.interiorEdges:
				a_CX[f02+e.indexAsInteriorEdge()] = a_lambda[e.index()];
			Coh1_CX.append(a_CX);
		Coh1_CX_matrix = matrix(Coh1_CX);
		self.Coh1_CX_matrix = Coh1_CX_matrix;
		self.Coh1_CX = Coh1_CX;

		#print 3;

		#Compute H_1(surface) wrt. CX, the edgeFlag-graph:
		C_CX = self.chainComplex_CX();
		H1_CX = [c[1].vector(1) for c in C_CX.homology(generators=True,deg=1)];
		H1_CX_matrix = matrix(H1_CX);
		#H1_CX = matrix(H1_CX).LLL().rows(); ###LLL not necessary because of the following step:
		#Choose basis that makes it dual to the cohomology basis:
		Trf = H1_CX_matrix * Coh1_CX_matrix.transpose();
		H1_CX_matrix = Trf.inverse() * H1_CX_matrix;
		H1_CX = H1_CX_matrix.rows();
		self.H1_CX_matrix = H1_CX_matrix;
		self.H1_CX = H1_CX;

		#print 4;

		#Compute H_1(surface) wrt. X, the vertexFlag-graph:
		H1_X = [];
		for c in H1_CX:
			c_X = vector([c[vf.index()] for vf in self.vertexFlags]);
			H1_X.append(c_X);
		H1_X_matrix = matrix(H1_X);
		self.H1_X_matrix = H1_X_matrix;
		self.H1_X = H1_X;

		#print 5;

		#Compute H^1(surface) wrt. X, the vertexFlag-graph:
		Coh1_X = [];
		for a_lambda in Coh1_lambda:
			a_X = zero_vector(f02);
			for vf in self.vertexFlags:
				if vf.isWhite():
					a_X[vf.index()] = a_lambda[vf.efLeft.edge.index()] - a_lambda[vf.efRight.edge.index()];
			Coh1_X.append(a_X);
		Coh1_X_matrix = matrix(Coh1_X);
		self.Coh1_X_matrix = Coh1_X_matrix;
		self.Coh1_X = Coh1_X;

	def deletionOfQuad(self,q):
		'''
		Returns a NEW discrete Riemann surface with quad q deleted.
		'''
		
		numQuads = self.f2()-1;
		iq = q.index();
		forbiddenEdgeIndices = [4*iq,4*iq+1,4*iq+2,4*iq+3];
		pairsOfIndicesOfGluedEdgeFlags = [];
		for e in self.edges:
			if e.isInteriorEdge():
				i1 = e.rightEdgeFlag().index();
				i2 = e.leftEdgeFlag().index();
				if i1>=4*iq:
					if i1 in forbiddenEdgeIndices:
						continue;
					i1 -= 4;
				if i2>=4*iq:
					if i2 in forbiddenEdgeIndices:
						continue;
					i2 -= 4;
				pairsOfIndicesOfGluedEdgeFlags.append((i1,i2));
		return DiscreteRiemannSurface(numQuads,pairsOfIndicesOfGluedEdgeFlags);
	
	def sumWithSurfaceAlongQuad(self,P2,q1=None,q2=None):
		'''
		Delete a face from each of the surfaces and then attach them:
		Important: q1 and q2 need to be INTERIOR faces of self and P2, respectively.
		'''
		
		if q1 == None:
			q1 = self.quads[-1];
		if q2 == None:
			q2 = P2.quads[-1];
		i0,i1,i2,i3 = [ef.oppositeEdgeFlag().index() for ef in q1.edgeFlags()];
		j0,j1,j2,j3 = [ef.oppositeEdgeFlag().index() for ef in q2.edgeFlags()];
		result = self.deletionOfQuad(q1);
		P2del = P2.deletionOfQuad(q2);
		result.attach(P2del,[(i0,j3),(i1,j2),(i2,j1),(i3,j0)]);
		return result;

	def slicingAlongEdge(self,e=None):
		'''
		Returns a NEW discrete Riemann surface obtained by slicing self along edge e.
		'''
		
		if e == None:
			e = self.interiorEdges[-1];
		pairsOfIndicesOfGluedEdgeFlags = [];
		for e1 in self.edges:
			if e1.isInteriorEdge() and e1 != e:
				i1 = e1.rightEdgeFlag().index();
				i2 = e1.leftEdgeFlag().index();
				pairsOfIndicesOfGluedEdgeFlags.append((i1,i2));
		return DiscreteRiemannSurface(self.f2(),pairsOfIndicesOfGluedEdgeFlags);

	def sumWithSurfaceAlongEdge(self,P2,e1=None,e2=None):
		'''
		Cut each the two surfaces along an edge and then glue them:
		Important: e1 and e2 need to be INTERIOR edges of self and P2, respectively.
		'''
		
		if e1 == None:
			e1 = self.interiorEdges[-1];
		if e2 == None:
			e2 = P2.interiorEdges[-1];
		i0,i1 = [ef.index() for ef in e1.edgeFlags];
		j0,j1 = [ef.index() for ef in e2.edgeFlags];
		result = self.slicingAlongEdge(e1);
		P2sliced = P2.slicingAlongEdge(e2);
		result.attach(P2sliced,[(i0,j1),(i1,j0)]);
		return result;

	def flippingAtVertex(self,v):
		if not v.isInteriorVertex():
			raise Exception("v needs to be an interior vertex!");
		efs = [];
		es = [];
		qs = [];
		for vf in v.vertexFlags:
			qs.append(vf.quad);
			e1 = vf.efRight.efRight();
			if not e1.edge.isInteriorEdge():
				raiseException("All edges in link of v need to be interior edges!");
			e2 = e1.efRight();
			if not e2.edge.isInteriorEdge():
				raiseException("All edges in link of v need to be interior edges!");
			efs.append(e1.oppositeEdgeFlag());
			efs.append(e2.oppositeEdgeFlag());
			es.append(e1.edge);
			es.append(e2.edge);
		if v.isWhite():
			ei = [3,0,1,2];
		else:
			ei = [0,1,2,3];
		gluings = [];
		#Remember all edges in the deletion of v:
		for e in self.interiorEdges:
			if not v in e.vertices():
				if e not in es:
					gluings.append((e.edgeFlags[0].index(),e.edgeFlags[1].index()));
		n = len(qs);
		for i in range(n):
			#Glue new quads to a cone around v:
			gluings.append((4*qs[i].index()+ei[3],4*qs[(i+1)%n].index()+ei[0]));
			#Glue new quads to old link of v:
			gluings.append((4*qs[i].index()+ei[1],efs[(2*i+1)%(2*n)].index()));
			gluings.append((4*qs[i].index()+ei[2],efs[(2*i+2)%(2*n)].index()));
		return DiscreteRiemannSurface(self.f2(),gluings);			

	def rhoStar(self,rho):
		return 1/rho;

	def qStar(self,q):
		return -1/q;

	def rhoOfQ(self,q):
		return I*(q+1)/(q-1);

	def qOfRho(self,rho):
		return (I*rho-1)/(I*rho+1);

	def trivialComplexStructure(self):
		'''
		Return the combinatorial complex structure where every
		quadrilateral is considered as a square.
		'''
		
		QLog = [(1/2,0) for i in range(self.f2())];
		Q = [exp(pi*I*theta + t) for theta,t in QLog];
		return QLog, Q;		

	def complexParallelogramStructure(self,theta0=[],t0=[],largeAngles=None):
		'''
		- theta0 is a list of (quad,theta), where quad is a quad of the mesh.
		  This prescribes the complex structure at quad.
		- t0 is a list of (quad,t), where quad is a quad of the mesh.
		  This prescribes the complex structure at quad.
		- largeAngles is a list of (v,k), where v is a vertex of the mesh.
		  This prescribes the conical angle at v as k*pi.
		  The other vertices will have conical angle 2*pi.
		  If largeAngles==None, then we will choose the conical angles
		    automatically w.r.t. a heuristic (and 2*pi mod 4*pi).
				
		Note that pi <-> 1, i.e.
		q = exp(pi*I*theta + t)
		
		We assume that surface is closed and has genus g>=1.
		'''
		
		#Equations around vertices:
		#v white: prod_quad sqrt(q) = -1:
		#			thus sum_quad 1/2*theta = angle
		#			and sum_quad t = 0	
		#v black: prod_quad sqrt(q^*) = -1
		#			thus sum_quad 1/2*(1-theta) = angle
		#			and sum_quad t = 0	

		angles = [2 for i in range(self.f0())];
		if largeAngles != None:
			for v,k in largeAngles:
				angles[v.index()] = k;
		else:
			remaining4pis = self.genus()-1; #Total curvature is 2pi times chi (and chi=2-2g).
			#Sort vertices of surface decreasingly by number of incident quads:
			vis = range(self.f0());
			while remaining4pis >= 1:
				vis.sort(key = lambda i: angles[i]/2-len(self.vertices[i].vertexFlags));
				angles[vis[0]] += 4;
				remaining4pis -= 1;	

		print "angles:",angles;
		
		eqnsTheta = [];
		rhsTheta = [];
		#Equations for prescribed thetas:
		for quad,theta in theta0:
			eq = zero_vector(QQ,self.f2());
			eq[quad.index()] = 1;
			rhs = theta;
			eqnsTheta.append(eq);
			rhsTheta.append(rhs);
		#Equations for thetas around vertices:
		for v in self.vertices:
			if v.isInteriorVertex():
				eq = zero_vector(QQ,self.f2());
				rhs = angles[v.index()];
				if v.isWhite():
					for quad in v.adjacentQuads():
						eq[quad.index()] += 1;
				else:
					for quad in v.adjacentQuads():
						eq[quad.index()] += -1;
						rhs -= 1;
				eqnsTheta.append(eq);
				rhsTheta.append(rhs);

		A = matrix(QQ,eqnsTheta);
		b = vector(QQ,rhsTheta);
		print A;
		print b;
		#print "For thetas: rank(A) =",A.rank()," and dim(b) =",len(b);
		Theta = A.solve_right(b); #Just to check solubility, TODO: rather check later A*Theta=b...
		#print "One possible not very good theta:",Theta;

		Astar = A.pseudoinverse();
		#print Astar;
		Theta = Astar * b; #Computes the Theta with the smallest 2-norm.
		#print "More balanced Theta:",Theta;

		eqnsT = [];
		rhsT = [];
		#Equations for prescribed ts:
		for quad,t in t0:
			eq = zero_vector(self.f2());
			eq[quad.index()] += 1;
			rhs = 0;
			eqnsT.append(eq);
			rhsT.append(rhs);
		#Equations for ts around vertices:
		for v in self.vertices:
			if v.isInteriorVertex():
				eq = zero_vector(self.f2());
				rhs = 0;
				for quad in v.adjacentQuads():
					eq[quad.index()] += 1;
				eqnsT.append(eq);
				rhsT.append(rhs);

		A = matrix(QQ,eqnsT);
		b = vector(QQ,rhsT);
		#print "For ts: rank(A) =",A.rank()," and dim(b) =",len(b);
		T = A.solve_right(b); #Just to check solubility, TODO: rather check later A*T=b...
		#print "T:",T;

		Astar = A.pseudoinverse();
		T = Astar * b; #Computes the T with the smallest 2-norm.

		QLog = [(Theta[i],T[i]) for i in range(self.f2())];
		Q = [exp(pi*I*theta + t) for theta,t in QLog];

		return QLog,Q;

	def randomComplexStructure(self):
		'''
		Note that pi <-> 1, i.e.
		q = exp(pi*I*theta + t)
		'''
		
		#Equations around vertices:
		#None
		
		QLog = [(random(),random()) for i in range(self.f2())];
		Q = [exp(pi*I*theta + t) for theta,t in QLog];

		return QLog,Q;

	def show(self,Q=None,quad0=None):
		'''
		Q is a complex Structure.
		'''
		
		if Q == None:
			QLog,Q = self.complexParallelogramStructure();
		if quad0 == None:
			quad0 = self.quads[0];
		#Breadth first seach:
		xyr = {quad0.index(): (CC(0),CC(1))};
		nextQuads = [quad0];
		usedEdges = [];
		while nextQuads != []:
			quad = nextQuads[0];
			nextQuads.remove(quad);
			q0 = Q[quad.index()]
			#Go through all adjacent quads:
			efs = quad.edgeFlags();
			for i in range(4):
				ef = efs[i];
				efOpp = ef.oppositeEdgeFlag();
				if efOpp != None:
					quadOpp = efOpp.quad;
					if not xyr.has_key(quadOpp.index()):
						qOpp = Q[quadOpp.index()];
						xy,r = xyr[quad.index()];
						toRotate = False;
						if i==0:
							xyr[quadOpp.index()] = (xy,r/qOpp);
							if quadOpp.ef3 != efOpp:
								toRotate = True;						
						elif i==1:
							xyr[quadOpp.index()] = (xy+r*(1+q0),-r*q0);
							if quadOpp.ef0 != efOpp:
								toRotate = True;						
						elif i==2:
							xyr[quadOpp.index()] = (xy+r*(1+q0),-r/qOpp);
							if quadOpp.ef3 != efOpp:
								toRotate = True;						
						elif i==3:
							xyr[quadOpp.index()] = (xy,r*q0);
							if quadOpp.ef0 != efOpp:
								toRotate = True;
						if toRotate:
							xy,r = xyr[quadOpp.index()];
							xyr[quadOpp.index()] = (xy+r*(1+qOpp),-r);
						nextQuads.append(quadOpp);
						usedEdges.append(ef.edge);
		G = Graphics();
		for quadi,(xy,r) in xyr.iteritems():
			q = Q[quadi];
			Ps = [CC(xy),CC(xy+r),CC(xy+r*(1+q)),CC(xy+r*q)];
			#G += polygon(Ps);
			G += line(Ps+[Ps[0]],zorder=-1);
			G += text(str(quadi),CC(xy+1/2*r*(1+q)),fontsize=14);


		for quadi,(xy,r) in xyr.iteritems():
			q = Q[quadi];
			quad = self.quads[quadi];
			Ps = [CC(xy),CC(xy+r),CC(xy+r*(1+q)),CC(xy+r*q)];
			for i in range(4):
				p = Ps[i];
				if i%2 == 0: #white:
					G += point(p,pointsize=100,rgbcolor=(0,0,0));
					G += point(p,pointsize=80,rgbcolor=(1,1,1));
				else: #black:
					G += point(p,pointsize=100,rgbcolor=(1,1,1));
					G += point(p,pointsize=80,rgbcolor=(0,0,0));
				if True: #quad.edge(i) not in usedEdges:
					xy = CC((Ps[i]+Ps[(i+1)%4])/2);
					G += point(xy,pointsize=140,rgbcolor=(1,1,1));
					G += text(str(quad.edge(i).index()),xy,fontsize=10);
		
		G.show(axes=False,aspect_ratio=1);

	def aCombinatorialSpinStructure(self):
		'''
		Computes a function mu from the vertexFlags to Z/2={+1,-1},
		such that integrating over local trivial cycles in
		the middle-edge-graph X yields -1.
		Such a function is called a combinatorial spin structure mu.
		In fact, adding to mu a coboundary in C^1(X;Z/2) yields an
		equivalent combinatorial spin structure.
		Moreover, H^1(X;Z/2) acts faithful and transitively on
		the set of spin structures.
		'''
		
		mu = [+1 for ef in self.edgeFlags];

		#First adjust mu such that it integrates to -1 around vertices:
		for v in self.vertices:
			if v.isInteriorVertex():
				mu[v.vertexFlags[0].efRight.index()] = -1;
		#print "After adjusting at vertices, mu:",mu;

		#Now adjust mu such that it integrates to -1 around quads:
		badQuads = set([]);
		for q in self.quads:
			muq = mu[q.ef0.index()] * mu[q.ef1.index()] * mu[q.ef2.index()] * mu[q.ef3.index()];
			if muq == +1:
				badQuads.add(q.index());
		#print "badQuads:",badQuads;
		Gdual = self.dualGraph_lambda();
		dist,pred = Gdual.shortest_path_all_pairs();
		while len(badQuads)>0:
			#Take two bad quads: #This works as len(badQuads) is even.
			q1 = self.quads[badQuads.pop()];
			q2 = self.quads[badQuads.pop()];
			while q1 != q2:
				q2pred = pred[q1][q2];
				#print q2pred,q2;
				e = q2.commonEdge(q2pred);
				#The following flips the badness of q2,
				#and also the badness of q2pred.
				#It does not change the badness around any vertex.
				mu[e.edgeFlags[0].vfWhite.index()] *= -1;
				mu[e.edgeFlags[1].vfWhite.index()] *= -1;
				q2 = q2pred;
		return mu;		

	def windingNumberAroundSpinStructure_X(self,mu_X,loop_X):
		'''
		A spin structure mu_X: vertexFlags -> {+1,-1} yields
		a tangent vector field with singularities of even index.
		1-cycles in X have thus a winding number mod 2 (in {0,1}) wrt. it.
		That's what this function is computing.
		Here, loop_X is a list of edgeFlags, which should represent a loop in X.
		'''
		
		windingNumber = 0; #we compute in Q/(2Z).
		for i in range(len(loop_X)):
			vf = loop_X[i];
			vf2 = loop_X[(i+1) % len(loop_X)];
			if mu_X[vf.index()] == -1:
				windingNumber += 1;
			if vf.quad == vf2.quad:
				if vf2 == vf.efLeft.theOtherVertexFlag(vf):
					windingNumber += 1/2;
				elif vf2 == vf.efRight.theOtherVertexFlag(vf):
					windingNumber -= 1/2;
				else:
					raise Exception("Given loop_X is not a loop!");		
		return Mod(ZZ(windingNumber),2);

	#Worked only for regular surfaces:
	def DEPRECIATED_johnsonFormAtCycle_X(self,mu_X,c_X):
		cMod2 = [Mod(ZZ(cvf),2) for cvf in c_X];

		nodes = copy(self.edges);
		edges = [(vf.efRight.edge,vf.efLeft.edge) for vf in self.vertexFlags if cMod2[vf.index()]!=0];
		cGraph = Graph([nodes,edges]);
		#cGraph.show();
		loopsInGraph = cGraph.cycle_basis(output="vertex");
		loops = [];
		loopsAsCycles = [];
		for l in loopsInGraph:
			c = [];
			cAsCycle = zero_vector(self.f02());
			for i in range(len(l)):
				e1 = l[i];
				e2 = l[(i+1) % len(l)];
				vf = e1.vertexFlagTowardsOtherEdge(e2);
				c.append(vf);
				cAsCycle[vf.index()] = +1;
			loops.append(c);
			loopsAsCycles.append(cAsCycle);
		
		result = Mod(0,2); #we compute mod 2:
		for c in loops:
			result += self.windingNumberAroundSpinStructure_X(mu_X,c) + 1;
		for i in range(len(loops)):
			for j in range(i+1,len(loops)):
				#Here it is important, that the loops are simple:
				result += self.intersectionNumber_X(loopsAsCycles[i],loopsAsCycles[j]);
		return result;

	def johnsonFormAtCycle_X(self,mu_X,c_X):
		'''
		Should work now also for non-regular surfaces!
		'''
		
		cMod2 = [Mod(ZZ(cvf),2) for cvf in c_X];
		remaining_vfs = [vf for vf in self.vertexFlags if cMod2[vf.index()] != 0];
		loops = [];
		loopsAsCycles = [];
		while len(remaining_vfs) > 0:
			c = [];
			cAsCycle = zero_vector(self.f02());
			vf = remaining_vfs[0];
			e0 = vf.efLeft.edge;
			ef = vf.efLeft;
			while True:
				remaining_vfs.remove(vf);
				c.append(vf);
				cAsCycle[vf.index()] += 1;
				ef = vf.theOtherEdgeFlag(ef);
				e = ef.edge;
				if e == e0: #closed the loop
					break;

				#Need to find the next vf in the loop:
				#For this we first check the adjacent vf in the same quad.
				#This ensures that the decomposition of cMod2 in loops
				#yields simple and almost disjoint loops
				#(disjoint after a small deformation).
				
				foundNewVertexFlag = False;
				efs = [ef];
				if e.isInteriorEdge():
					efs.append(ef.oppositeEdgeFlag());
				for efNew in efs: 
					for vfNew in efNew.vertexFlags():
						if vfNew in remaining_vfs:
							vf = vfNew;
							foundNewVertexFlag = True;
							break;
					if foundNewVertexFlag:
						ef = efNew;
						break;
				if not foundNewVertexFlag:
					raise Exception("Given chain c_X is mod-2 not a cycle!");
			loops.append(c);
			loopsAsCycles.append(cAsCycle);
		
		result = Mod(0,2); #we compute mod 2:
		for c in loops:
			result += self.windingNumberAroundSpinStructure_X(mu_X,c) + 1;
		for i in range(len(loops)):
			for j in range(i+1,len(loops)):
				#Here it is important, that the loops are simple:
				result += self.intersectionNumber_X(loopsAsCycles[i],loopsAsCycles[j]);
		return result;

	def arfInvariant(self,mu_X=None):
		'''
		Compute the Arf invariant of the Johnson quadratic form
		associated to a Spin structure mu_X: edgeFlags -> {-1,+1}.
		'''
		
		if mu_X == None:
			mu_X = self.aCombinatorialSpinStructure();
		if self.H1_X == None:
			self.computeHomologies();
		g = self.genus();
		twoG = len(self.H1_X); #= 2*g
		sumQ = 0;
		for w in Words([0,1],twoG):
			#print '====================================';
			lamda = vector(w);
			#print lamda;
			gamma_X = lamda * self.H1_X_matrix;
			#print gamma_X;
			qLamda = self.johnsonFormAtCycle_X(mu_X,gamma_X);
			#qLamda = self.DEPRECIATED_johnsonFormAtCycle_X(mu_X,gamma_X);
			#print qLamda;
			
			sumQ += qLamda.lift(); #lfit makes it to an element 0 or 1 in Z.
		#print sumQ;
		if sumQ == 2^(2*g-1) + 2^(g-1):
			return 1;
		elif sumQ == 2^(2*g-1) - 2^(g-1):
			return 0;
		else:
			raise Exception("Something is wrong when computing the Arf invariant!");

	def translateCombinatorialSpinStructureViaCohomologyClass(self,mu_X,alpha_X):
		'''
		Here, mu_X: edgeFlags -> {+1,-1} represents a spin structure, and
		alpha_X: edgeFlags -> {0,1} represents an element in H^1(X;F_2).
		'''
		
		mu_X_new = [ZZ(mu_X[ef.index()] * (-1)^abs(alpha_X[ef.index()])) for ef in self.edgeFlags];
		return mu_X_new;

	def spinStructures(self):
		mu = self.aCombinatorialSpinStructure();
		if self.Coh1_X_matrix == None:
			self.computeHomologies();
		g = self.genus();
		Mus = [];
		for w in Words([0,1],2*g):
			coeff = vector(w);
			alpha = coeff * self.Coh1_X_matrix;
			mu_coeff = self.translateCombinatorialSpinStructureViaCohomologyClass(mu,alpha);
			Mus.append(mu_coeff);
		return Mus;

	def testArfInvariants(self,q=None):
		'''
		q: quads -> C is a combinatorial complex structure.
		'''
		
		if q == None:
			qLog, q = self.trivialComplexStructure();
			#qLog, q = self.complexParallelogramStructure();
			#qLog, q = self.complexParallelogramStructure(theta0=[(self.quads[0],1/5)]);
			#qLog, q = self.randomComplexStructure();
		mu = self.aCombinatorialSpinStructure();
		if self.Coh1_X_matrix == None:
			self.computeHomologies();
		g = self.genus();
		sumArfs = 0;
		for w in Words([0,1],2*g):
			coeff = vector(w);
			print "======== coeff =",coeff;
			alpha = coeff * self.Coh1_X_matrix;
			mu_coeff = self.translateCombinatorialSpinStructureViaCohomologyClass(mu,alpha);
			arf = self.arfInvariant(mu_coeff);
			print "arf =",arf;
			sumArfs += arf;
			M_CC = self.holomorphicVSpinorsMatrix(q=q,mu_X=mu_coeff,K=CC);
			M_CIF = self.holomorphicVSpinorsMatrix(q=q,mu_X=mu_coeff,K=CIF);
			print "-"
			if not M_CIF.det().contains_zero():
				print "nullity(M) = 0 (provably)";
			else:
				eigenvalues = M_CC.eigenvalues();
				eigenvalues_abs = [abs(lamda) for lamda in eigenvalues];
				eigenvalues_abs.sort();
				print "eigenvalues_abs(M):",eigenvalues_abs;
		print "Ratio of arf-invariant one:",sumArfs/2^(2*g);

	def holomorphicVSpinorsMatrix(self,q=None,mu_X=None,K=CIF):
		'''
		q: quads -> C is a combinatorial complex structure.
		mu_X: vertexFlags -> {+1,-1} is a combinatorial Spin structure.
		K is the base field for the computation.
		This computes the square matrix whose kernel corresponds to
		holomorphic V-spinors.
		'''
		
		if q == None:
			qLog, q = self.trivialComplexStructure();
			#qLog, q = self.complexParallelogramStructure();
			#qLog, q = self.complexParallelogramStructure(theta0=[(self.quads[0],1/5)]);
		if mu_X == None:
			mu_X = self.aCombinatorialSpinStructure();

		eqns = [];
		for quad in self.quads:
			rhoQuad = self.rhoOfQ(q[quad.index()]);
			rhoQuadStar = self.rhoStar(rhoQuad);
			#First equation:
			eqn = zero_vector(K,self.f1());
			eqn[quad.ef2.edge.index()] = mu_X[quad.vf2.index()];
			eqn[quad.ef1.edge.index()] = -sqrt(1+rhoQuad^2);
			eqn[quad.ef0.edge.index()] = mu_X[quad.vf1.index()] * rhoQuad;
			eqns.append(eqn);
			#Second equation:
			eqn = zero_vector(K,self.f1());
			eqn[quad.ef3.edge.index()] = mu_X[quad.vf3.index()];
			eqn[quad.ef2.edge.index()] = -sqrt(1+rhoQuadStar^2);
			eqn[quad.ef1.edge.index()] = mu_X[quad.vf2.index()] * rhoQuadStar;
			eqns.append(eqn);
			#There are in total 4 equations, but two of them being redundant.

		M = matrix(K,eqns);
		#print "det(M) =",det(M);
		#if K==CC:
		#	print "nullity(M) =",M.nullity(); #Has numerical issues!!!
		#if det(M).contains_zero():
		#	print M.nullity();
		#print "arf(mu) =",self.arfInvariant(mu_X);
		return M;

########################################################################
############# Particular surfaces: #####################################
########################################################################

def xOrigami(genus):
	'''
	Computes a non-regular surface of genus g with 2g quads.
    
	   ###2#####1#####4#####3################
	   #     #     #     #     #            #
	   A  0  #  1  #  2  #  3  #    ....    A
	   #     #     #     #     #            #
	   ###1#####2#####3#####4################
	'''
	
	#Former implementation for g=3:
	#return DiscreteRiemannSurface(6,[(1,6),(4,11),(9,14),(12,19),(17,22),(3,20),(0,5),(2,7),(8,13),(10,15),(16,21),(18,23)]);

	f2 = genus*2;
	gluings = [];
	for i in range(genus):
		q0 = 2*i;
		gluings.append((4*q0+0,4*q0+5));
		gluings.append((4*q0+2,4*q0+7));
		gluings.append((4*q0+1,4*q0+6));
		gluings.append((4*q0+4,(4*q0+3+8) % (f2*4)));
	return DiscreteRiemannSurface(f2,gluings);

def xOrigami2(genus):
	'''
	Computes an almost regular surface of genus g with 4g quads.
    
	   #########2#####1#####4#####3##############
	   #     #     #     #     #     #          #
	   A  2  #  3  #  6  #  7  #     #   ....   A
	   #     #     #     #     #     #          #
	   ##########################################
	   #     #     #     #     #     #          #
	   B  0  #  1  #  4  #  5  #     #   ....   B
	   #     #     #     #     #     #          #
	   ###1#####2#####3#####4####################
	'''
	
	f2 = genus*4;
	gluings = [];
	for i in range(genus):
		q0 = 4*i;
		gluings.append((4*q0+0,(4*q0+9+16) % (f2*4)));
		gluings.append((4*q0+7,4*q0+14));
		gluings.append((4*q0+1,4*q0+6));
		gluings.append((4*q0+5,4*q0+12));
		gluings.append((4*q0+8,4*q0+15));
		gluings.append((4*q0+2,4*q0+11));
		gluings.append((4*q0+4,(4*q0+3+16) % (f2*4)));
		gluings.append((4*q0+13,(4*q0+10+16) % (f2*4)));
	return DiscreteRiemannSurface(f2,gluings);

def xOrigami3(genus):
	'''
	Computes another almost regular surface of genus g with 6g quads.
	
	   ###2#####1#####4#####3##############
	   #     #     #     #     #          #
	   A  4  #  5  # 10  # 11  #   ....   A
	   #     #     #     #     #          #
	   ####################################
	   #     #     #     #     #          #
	   B  2  #  3  #  8  #  9  #   ....   B
	   #     #     #     #     #          #
	   ####################################
	   #     #     #     #     #          #
	   C  0  #  1  #  6  #  7  #   ....   C
	   #     #     #     #     #          #
	   ###1#####2#####3#####4##############
	'''
	
	f2 = genus*6;
	gluings = [];
	for i in range(genus):
		q0 = 6*i;
		gluings.append((4*q0+0,4*q0+21));
		gluings.append((4*q0+7,4*q0+18));
		gluings.append((4*q0+1,4*q0+6));
		gluings.append((4*q0+5,4*q0+12));
		gluings.append((4*q0+8,4*q0+15));
		gluings.append((4*q0+2,4*q0+11));
		gluings.append((4*q0+9,4*q0+16));
		gluings.append((4*q0+14,4*q0+23));
		gluings.append((4*q0+17,4*q0+22));
		gluings.append((4*q0+4,(4*q0+3+24) % (f2*4)));
		gluings.append((4*q0+13,(4*q0+10+24) % (f2*4)));
		gluings.append((4*q0+20,(4*q0+19+24) % (f2*4)));
	return DiscreteRiemannSurface(f2,gluings);

def xOrigami4(genus):
	'''
	Computes a regular surface of genus g with 12g quads.
	It is regular, and has less faces than the similar
	regular genus g surface xOrigami(g).subdivision(3).
	
	   #########3#####4#####1#####2########
	   #     #     #     #     #          #
	   A  4  #  5  # 10  # 11  #   ....   A
	   #     #     #     #     #          #
	   ####################################
	   #     #     #     #     #          #
	   B  2  #  3  #  8  #  9  #   ....   B
	   #     #     #     #     #          #
	   ####################################
	   #     #     #     #     #          #
	   C  0  #  1  #  6  #  7  #   ....   C
	   #     #     #     #     #          #
	   ###1#####2#####3#####4##############
	'''
	
	f2 = genus*12;
	gluings = [];
	for i in range(2*genus):
		q0 = 6*i;
		if i%2 == 0:
			gluings.append((4*q0+0,(4*q0+21+24) % (f2*4)));
			gluings.append((4*q0+7,(4*q0+18+48) % (f2*4)));
		else:
			gluings.append((4*q0+0,(4*q0+21-24) % (f2*4)));
			gluings.append((4*q0+7,(4*q0+18+0) % (f2*4)));
		gluings.append((4*q0+1,4*q0+6));
		gluings.append((4*q0+5,4*q0+12));
		gluings.append((4*q0+8,4*q0+15));
		gluings.append((4*q0+2,4*q0+11));
		gluings.append((4*q0+9,4*q0+16));
		gluings.append((4*q0+14,4*q0+23));
		gluings.append((4*q0+17,4*q0+22));
		gluings.append((4*q0+4,(4*q0+3+24) % (f2*4)));
		gluings.append((4*q0+13,(4*q0+10+24) % (f2*4)));
		gluings.append((4*q0+20,(4*q0+19+24) % (f2*4)));
	return DiscreteRiemannSurface(f2,gluings);

def xOrigami5(genus):
	'''
	Computes a regular surface of genus g with 8g quads.
	It is regular, and has less faces than the similar
	regular genus g surface xOrigami(g).subdivision(3).
	
	   ###3#####4#####1#####2##############
	   #     #     #     #     #          #
	   A  2  #  3  #  7  #  5  #   ....   A
	   #     #     #     #     #          #
	   ####################################
	   #     #     #     #     #          #
	   C  0  #  1  #  6  #  4  #   ....   C
	   #     #     #     #     #          #
	   ###1#####2#####3#####4##############
	'''
	
	return xOrigami(genus).subdivision(2);

def quaternionOrigami():
	'''
	Create a surface with one quad per element in the quaternion group Q.
	For g in Q, quad g is glued to the left of g*i, and below g*j.
	To match vertex colors, quads are rotated by 90 degrees if
	 the corresponding element g maps to the non-zero element in Q/<k>.
	
	         ...
	
	       #######
	       #     #
	   ... # g*j # ...
	       #     #
	       #############
	       #     #     #
	   ... #  g  # g*i # ...
	       #     #     #
	       #############
	
	         ...   ...
	
	  i.e.
	
	                     ###2#####1#####4#####3###
	                     #     #     #     #     #
	                     9 -k  # -j  #  k  #  j  9
	                     #     #     #     #     #
	   ###7#####6#####5###########5#####6#####7###
	   #     #     #     #     #
	   8  1  #  i  # -1  # -i  8
	   #     #     #     #     #
	   ###1#####2#####3#####4###
	'''
	
	Q = QuaternionGroup();
	L = list(Q);
	i = Q.gen(0);
	j = Q.gen(1);
	k = i*j;
	S = [g for g in Q if g in Q.subgroup([k])];
	ei = {True: [0,1,2,3], False: [3,0,1,2]}
	gluings = [];
	for g in Q:
		gluings.append((4*L.index(g)+ei[g in S][1],4*L.index(g*i)+ei[g*i in S][3]));
		gluings.append((4*L.index(g)+ei[g in S][2],4*L.index(g*j)+ei[g*j in S][0]));
	#print gluings;
	return DiscreteRiemannSurface(8,gluings);	

#Doesn't work as graph is not bipartite:
def DEPRECIATED_starOrigami(genus):
	'''
	Creates a star origami, a closed surface of genus g with 2g quads.
	Here, genus needs to be even to be indeed the genus!
	'''
	
	gluings = [];
	f2 = 2*genus;
	for q in range(f2-1):
		e0 = 4*q;
		if q % 2 == 0:
			gluings.append((e0+1,e0+2+4));
			gluings.append((e0+3,e0+0+4));
		else:
			gluings.append((e0+1,e0+0+4));
			gluings.append((e0+3,e0+2+4));
	gluings.append((0,2));
	gluings.append((4*f2-3,4*f2-1));
	return DiscreteRiemannSurface(f2,gluings);

def starOrigami_MoreSymmetric(oddGenus):
	'''
	Creates a version of the star origami, a closed surface
	of genus g (if g is ODD) with 2g quads.
	Here, genus needs to be even to be indeed the genus!
	
	             ###3######
	             #     #
	             #2i+2 # ...
	             #     #
	       ###1############
	       #     #     #
	       2 2i  #2i+1 2 
	       #     #     #
	    ############3###
	       #     #
	   ... #2i-1 #
	       #     #
	    ######1###
	'''
	
	gluings = [];
	g = oddGenus;
	f2 = 2*g;
	for q in range(f2):
		e0 = 4*q;
		if q % 2 == 0:
			gluings.append((e0+1,(e0+2+4) % (4*f2)));
			gluings.append((e0+3,(e0+0+4) % (4*f2)));
		else:
			gluings.append((e0+1,(e0+0+4) % (4*f2)));
			gluings.append((e0+3,(e0+2+4) % (4*f2)));
	return DiscreteRiemannSurface(f2,gluings);

def halfCube():
	return DiscreteRiemannSurface(3,[(3,4),(7,8),(0,11)]);

def cube():
	return DiscreteRiemannSurface(6,[(3,4),(7,8),(0,11),(2,15),(5,14),(6,19),(9,18),(10,23),(1,22),(13,16),(17,20),(21,12)]);

def cubeWithoutNFacets(n):
	C = cube();
	for i in range(n):
		C = C.deletionOfQuad(q=C.quads[-1]);
	return C;

def torus(n,m):
	'''
	Torus with 2n x 2m quads. Need n,m >= 1.
	
	   ###1#############n###
	   #     #       #     #
	  A1     #  ...  #     A1
	   #     #       #     #
	   #####################
	   #     #       #     #
	   #  :  #   :   #  :  #
	   #  :  #  ...  #  :  #
	   #  :  #   :   #  :  #
	   #     #       #     #
	   #####################
	   #     #       #     #
	  An     #  ...  #     An
	   #     #       #     #
	   ###1#############n###
	'''
	
	f1 = 4 * 2*n * 2*m;
	f2 = 2*n * 2*m;
	gluing = [];
	for i in range(2*n):
		for j in range(2*m):
			if (i+j)%2 == 0: 
				q = i + 2*n*j;
				if i<2*n-1:
					q_right = q + 1;
				else:
					q_right = q - (2*n-1);
				if i>0:
					q_left = q - 1;
				else:
					q_left = q + (2*n-1);
				if j<2*m-1:
					q_top = q + 2*n;
				else:
					q_top = q - 2*n*(2*m-1);
				if j>0:
					q_bottom = q - 2*n;
				else:
					q_bottom = q + 2*n*(2*m-1);
				gluing.append((4*q+0,4*q_bottom+1));
				gluing.append((4*q+1,4*q_right +2));
				gluing.append((4*q+2,4*q_top   +3));
				gluing.append((4*q+3,4*q_left  +0));
	#print gluing;
	return DiscreteRiemannSurface(f2,pairsOfIndicesOfGluedEdgeFlags=gluing);

def genus_n_surface(n,torusSize=1):
	'''
	Gluing n tori minus one quad.
	'''
	
	if n==0:
		return cube();
	elif n==1:
		return torus(torusSize,torusSize);
	elif n>=2:
		return genus_n_surface(n-1).sumWithSurfaceAlongQuad(torus(torusSize,torusSize));
	else:
		raise Exception();

def genus_n_surface_new_nonRegular(n,torusSize=1):
	'''
	Gluing n tori sliced along an edge.
	'''
	
	if n==0:
		return cube();
	elif n==1:
		return torus(torusSize,torusSize);
	elif n>=2:
		P_nm1 = genus_n_surface_new_nonRegular(n-1);
		T = torus(torusSize,torusSize);
		#Glue P_nm1 and T along the edge of P_nm1 with highest valence:
		e1 = max(P_nm1.edges,key = lambda e: len(e.adjacentEdges()));
		return P_nm1.sumWithSurfaceAlongEdge(T,e1=e1);
	else:
		raise Exception();

def genus_n_surface_new(n,torusSize=1):
	'''
	Glue n tori sliced along an edge, and then subdivide.
	'''
	
	return genus_n_surface_new_nonRegular(n,torusSize).subdivision(2);

########################################################################
############# Tests: ###################################################
########################################################################

#P = halfCube();
#P = cube();
#P = cubeWithoutNFacets(2);
#P = torus(1,1);
#P = torus(2,2);
#P = torus(3,3);
#P = genus_n_surface(2);
#P = genus_n_surface(2,torusSize=2);
#P = genus_n_surface(3);
#P = genus_n_surface(4);
#P = xOrigami(2);
#P = xOrigami(2).subdivision(4);
#P = xOrigami(3);
#P = xOrigami(3).subdivision(3);
#P = xOrigami(5).subdivision();
#P = xOrigami3(3)
#P = xOrigami(3).subdivision(3);
#P = genus_n_surface_new_nonRegular(2);
#P = genus_n_surface_new(3);
#P = quaternionOrigami();
#P = quaternionOrigami().subdivision(2);
#P = quaternionOrigami().subdivision(3);
#P = starOrigami(3);
#P = starOrigami_MoreSymmetric(3);

#P = P.flippingAtVertex(P.vertices[0]);
#print P;
#P.graph_lambda().show();

def testIntersectionForm(P):
	C = P.chainComplex_lambda();
	#print P.homology();
	print P.homologyBasis_lambda();
	print P.intersectionForm_lambda();
	print P.normalizedHomologyBasis_lambda();
	print P.intersectionForm(P.normalizedHomologyBasis_lambda());

#testIntersectionForm(P);

#print P.complexParallelogramStructure();
#print P.complexParallelogramStructure(theta0=[(P.quads[0],1/3)]);

#P.show();
#P.show(Q=P.complexParallelogramStructure(theta0=[(P.quads[0],1/5)])[1]);

#print P.aCombinatorialSpinStructure();

#hB = P.normalizedHomologyBasis_lambda();
#cohB = P.cohomologyBasis_lambda(hB);
#print cohB;

#P.computeHomologies();
#print P.arfInvariant();
#P.testArfInvariants();
#P.testArfInvariants(q=P.complexParallelogramStructure(theta0=[(P.quads[0],1/5)])[1]);
#P.vertexDegrees_edges();
