/******************************************************************
typePatchConsSparseSorted.h
*******************************************************************/

#ifndef __TYPEPatchConsSparseSorted_H__
#define __TYPEPatchConsSparseSorted_H__

#include <string.h>
#include <assert.h>
#include <limits> 

template <class T> class MRFEnergy;

class TypePartialEnumeration
{
private:
	struct Vector; // node parameters and messages
	struct Edge; // stores edge information and either forward or backward message

public:
	typedef enum
	{
		PatchConsSparseSortedVert,
		PatchConsSparseSortedHor,
	} Type;

	// types declarations
	typedef int Label;
	typedef double REAL;
	struct GlobalSize; // global information about number of labels
	struct LocalSize; // local information about number of labels (stored at each node)
	struct NodeData; // argument to MRFEnergy::AddNode()
	struct EdgeData; // argument to MRFEnergy::AddEdge()

	struct GlobalSize
	{
		GlobalSize(	int K, 
					int patchSize,
					int * labels,
					int * hor_source,
					int * hor_sink,
					int * vert_source,
					int * vert_sink,
					bool sort);

		int m_K;
		int m_patchSize;
		int * m_labels;
		int * m_hor_source;
		int * m_hor_sink;
		int * m_vert_source;
		int * m_vert_sink;
	};

	struct LocalSize // number of labels is stored at MRFEnergy::m_Kglobal
	{
		LocalSize(int K);
	};

	struct NodeData
	{
		NodeData(REAL* data); // data = pointer to array of size MRFEnergy::m_Kglobal

	private:
	friend struct Vector;
	friend struct Edge;
		REAL*		m_data;
	};

	struct EdgeData
	{
		EdgeData(Type type); // type must be PatchConsSparseSorted.

	private:
	friend struct Vector;
	friend struct Edge;
		Type		m_type;
	};

	static void SortVertical(int* labels, int* sortedLabels, int num_labels, int patch_size,bool issource);
	static void SortHorizontal(int* labels, int* sortedLabels, int num_labels, int patch_size,bool issource);

	static int twopow(int K)
	{
		int i, prod = 1;

		for (i=0;i<K;i++)
			prod = 2*prod;

		return prod;
	}

	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////// Visible only to MRFEnergy /////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

private:
friend class MRFEnergy<TypePartialEnumeration>;

	struct Vector
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize K); // returns -1 if invalid K's
		void Initialize(GlobalSize Kglobal, LocalSize K, NodeData data);  // called once when user adds a node
		void Add(GlobalSize Kglobal, LocalSize K, NodeData data); // called once when user calls MRFEnergy::AddNodeData()

		void SetZero(GlobalSize Kglobal, LocalSize K);                            // set this[k] = 0
		void Copy(GlobalSize Kglobal, LocalSize K, Vector* V);                    // set this[k] = V[k]
		void Add(GlobalSize Kglobal, LocalSize K, Vector* V);                     // set this[k] = this[k] + V[k]
		REAL GetValue(GlobalSize Kglobal, LocalSize K, Label k);                  // return this[k]
		REAL ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin);            // return min_k { this[k] }, set kMin
		REAL ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K);              // same as previous, but additionally set this[k] -= vMin (and kMin is not returned)

		static int GetArraySize(GlobalSize Kglobal, LocalSize K);
		REAL GetArrayValue(GlobalSize Kglobal, LocalSize K, int k); // note: k is an integer in [0..GetArraySize()-1].
		                                                            // For Potts functions GetArrayValue() and GetValue() are the same,
		                                                            // but they are different for, say, 2-dimensional labels.
		void SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x);

	private:
	friend struct Edge;
		REAL		m_data[1]; // actual size is MRFEnergy::m_Kglobal
	};

	struct Edge
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data); // returns -1 if invalid data
		static int GetBufSizeInBytes(int vectorMaxSizeInBytes); // returns size of buffer need for UpdateMessage()
		void Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj); // called once when user adds an edge
		Vector* GetMessagePtr();
		void Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj); // if the client calls this function, then the meaning of 'dir'
								                                               // in distance transform functions is swapped

		// When UpdateMessage() is called, edge contains message from dest to source.
		// The function must replace it with the message from source to dest.
		// The update rule is given below assuming that source corresponds to tail (i) and dest corresponds
		// to head (j) (which is the case if dir==0).
		//
		// 1. Compute Di[ki] = gamma*source[ki] - message[ki].  (Note: message = message from j to i).
		// 2. Compute distance transform: set
		//       message[kj] = min_{ki} (Di[ki] + V(ki,kj)). (Note: message = message from i to j).
		// 3. Compute vMin = min_{kj} m_message[kj].
		// 4. Set m_message[kj] -= vMin.
		// 5. Return vMin.
		//
		// If dir==1 then source corresponds to j, sink corresponds to i. Then the update rule is
		//
		// 1. Compute Dj[kj] = gamma*source[kj] - message[kj].  (Note: message = message from i to j).
		// 2. Compute distance transform: set
		//       message[ki] = min_{kj} (Dj[kj] + V(ki,kj)). (Note: message = message from j to i).
		// 3. Compute vMin = min_{ki} m_message[ki].
		// 4. Set m_message[ki] -= vMin.
		// 5. Return vMin.
		//
		// If Edge::Swap has been called odd number of times, then the meaning of dir is swapped.
		//
		// Vector 'source' must not be modified. Function may use 'buf' as a temporary storage.
		REAL UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* buf);


		// If dir==0, then sets dest[kj] += V(ksource,kj).
		// If dir==1, then sets dest[ki] += V(ki,ksource).
		// If Swap() has been called odd number of times, then the meaning of dir is swapped.
		void AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir);

	private:
		REAL horintcost(int i, int j, int dir, int K, int patchSize);
		REAL vertintcost(int i, int j, int dir, int K, int patchSize);

	protected:

		Type		m_type;

		// message
		Vector*		m_message;
	};

	struct EdgePatchConsSparseSorted : Edge
	{
	private:
	friend struct Edge;
		int		m_dir; // 0 if Swap() was called even number of times, 1 otherwise
	};
};


//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


inline TypePartialEnumeration::LocalSize::LocalSize(int K)
{
}

inline TypePartialEnumeration::GlobalSize::GlobalSize(	int num_labels, 
															int patchSize,
															int * labels,
															int * hor_source,
															int * hor_sink,
															int * vert_source,
															int * vert_sink,
															bool sort = true)
{
	m_K = num_labels;	
	m_patchSize = patchSize;
	m_labels = labels;

	// For each label set the sorting will be identical this can be done offline to speed up.
	if (sort)
	{
	 	SortHorizontal(labels, hor_source, num_labels, patchSize, true);
	 	SortHorizontal(labels, hor_sink, num_labels, patchSize, false);
	 	SortVertical(labels, vert_source, num_labels, patchSize, true);
	 	SortVertical(labels, vert_sink, num_labels, patchSize, false);
	}

	m_hor_source = hor_source;
	m_hor_sink = hor_sink;
	m_vert_source = vert_source;
	m_vert_sink = vert_sink;
}

/// Sorting the labels for O(n) message updates ///
inline void TypePartialEnumeration::SortHorizontal(int* labels, 
																								 			int* sortedLabels, 
																								 			int num_labels, 
																								 			int patchSize,
																								 			bool issource)
{
	int i, mask;
	int* lastPos = sortedLabels;
	int j;

	mask = twopow(patchSize*patchSize)-1;
	for (i=0;i<patchSize;i++)
		mask = mask & (~twopow((i*patchSize)));

	if (issource)
	{
		for (i = 0; i < twopow(patchSize*patchSize) ; i++)
		{
			if ((mask & i) == i) // Search for all patches consistent with this patch.
			{
				for (j=0; j<num_labels; j++)
				{
					if ((((int)labels[j]) & mask) == i)
					{
						*lastPos = j;
						++lastPos;
					}
				}
			}
		}
	}
	else
	{
		mask = mask >> 1;
		for (i = 0; i < twopow(patchSize*patchSize); i++)
		{
			if ((mask & i) == i) // Search for all patches consistent with this patch.
			{
				for (j=0; j<num_labels; j++)
				{
					if ((((int)labels[j]) & mask) == i)
					{
						*lastPos = j;
						++lastPos;
					}
				}
			}
		}
	}
}


inline void TypePartialEnumeration::SortVertical(int* labels, 
																						  int* sortedLabels, 
																						 	int num_labels, 
																						 	int patchSize,
																						 	bool issource)
{
	int i, mask;
	int* lastPos = sortedLabels;
	int j;

	mask = twopow(patchSize*patchSize)-1;
	for (i=0;i<patchSize;i++)
		mask = mask & (~twopow(i));


	if (issource)
	{
		for (i = 0; i < twopow(patchSize*patchSize); i++)
		{
			if ((mask & i) == i) // Search for all patches consistent with this patch.
			{
				for (j=0; j<num_labels; j++)
				{
					if ((((int)labels[j]) & mask) == i)
					{
						*lastPos = j;
						++lastPos;
					}
				}
			}
		}
	}
	else
	{
		mask = mask >> patchSize;
		for (i = 0; i < twopow(patchSize*patchSize); i++)
		{
			if ((mask & i) == i) // Search for all patches consistent with this patch.
			{
				for (j=0; j<num_labels; j++)
				{
					if ((((int)labels[j]) & mask) == i)
					{
						*lastPos = j;
						++lastPos;
					}
				}
			}
		}
	}
}

///////////////////// NodeData and EdgeData ///////////////////////

inline TypePartialEnumeration::NodeData::NodeData(REAL* data)
{
	m_data = data;
}


inline TypePartialEnumeration::EdgeData::EdgeData(Type type)
{
	assert((type == PatchConsSparseSortedHor) | (type == PatchConsSparseSortedVert));
	m_type = type;
}

///////////////////// Vector ///////////////////////

inline int TypePartialEnumeration::Vector::GetSizeInBytes(GlobalSize Kglobal, LocalSize K)
{
	if (Kglobal.m_K < 1)
	{
		return -1;
	}
	return Kglobal.m_K*sizeof(REAL);
}

inline void TypePartialEnumeration::Vector::Initialize(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	memcpy(m_data, data.m_data, Kglobal.m_K*sizeof(REAL));
}

inline void TypePartialEnumeration::Vector::Add(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	for (int k=0; k<Kglobal.m_K; k++)
	{
		m_data[k] += data.m_data[k];
	}
}

inline void TypePartialEnumeration::Vector::SetZero(GlobalSize Kglobal, LocalSize K)
{
	memset(m_data, 0, Kglobal.m_K*sizeof(REAL));
}

inline void TypePartialEnumeration::Vector::Copy(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	memcpy(m_data, V->m_data, Kglobal.m_K*sizeof(REAL));
}

inline void TypePartialEnumeration::Vector::Add(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	for (int k=0; k<Kglobal.m_K; k++)
	{
		m_data[k] += V->m_data[k];
	}
}

inline TypePartialEnumeration::REAL TypePartialEnumeration::Vector::GetValue(GlobalSize Kglobal, LocalSize K, Label k)
{
	assert(k>=0 && k<Kglobal.m_K);
	return m_data[k];
}

inline TypePartialEnumeration::REAL TypePartialEnumeration::Vector::ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin)
{
	REAL vMin = m_data[0];
	kMin = 0;
	for (int k=1; k<Kglobal.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
			kMin = k;
		}
	}

	return vMin;
}

inline TypePartialEnumeration::REAL TypePartialEnumeration::Vector::ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K)
{
	REAL vMin = m_data[0];
	for (int k=1; k<Kglobal.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
		}
	}
	for (int k=0; k<Kglobal.m_K; k++)
	{
		m_data[k] -= vMin;
	}

	return vMin;
}

inline int TypePartialEnumeration::Vector::GetArraySize(GlobalSize Kglobal, LocalSize K)
{
	return Kglobal.m_K;
}

inline TypePartialEnumeration::REAL TypePartialEnumeration::Vector::GetArrayValue(GlobalSize Kglobal, LocalSize K, int k)
{
	assert(k>=0 && k<Kglobal.m_K);
	return m_data[k];
}

inline void TypePartialEnumeration::Vector::SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x)
{
	assert(k>=0 && k<Kglobal.m_K);
	m_data[k] = x;
}

///////////////////// EdgeDataAndMessage implementation /////////////////////////

inline int TypePartialEnumeration::Edge::GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data)
{
	int messageSizeInBytes = Kglobal.m_K*sizeof(REAL);

	switch (data.m_type)
	{
		case PatchConsSparseSortedHor:
			return sizeof(EdgePatchConsSparseSorted) + messageSizeInBytes;
		case PatchConsSparseSortedVert:
			return sizeof(EdgePatchConsSparseSorted) + messageSizeInBytes;
		default:
			return -1;
	}
}

inline int TypePartialEnumeration::Edge::GetBufSizeInBytes(int vectorMaxSizeInBytes)
{
	return vectorMaxSizeInBytes;
}

inline void TypePartialEnumeration::Edge::Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj)
{
	m_type = data.m_type;
	int* dataptr;

	switch (m_type)
	{
		case PatchConsSparseSortedHor:
			((EdgePatchConsSparseSorted*)this)->m_dir = 0;
			m_message = (Vector*)((char*)this + sizeof(EdgePatchConsSparseSorted));
			break;
		case PatchConsSparseSortedVert:
			((EdgePatchConsSparseSorted*)this)->m_dir = 0;
			m_message = (Vector*)((char*)this + sizeof(EdgePatchConsSparseSorted));
			break;
		default:
			assert(0);
	}

	memset(m_message->m_data, 0, Kglobal.m_K*sizeof(REAL));
}


inline TypePartialEnumeration::Vector* TypePartialEnumeration::Edge::GetMessagePtr()
{
	return m_message;
}

inline void TypePartialEnumeration::Edge::Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj)
{
	if (m_type == PatchConsSparseSortedHor | m_type == PatchConsSparseSortedVert)
	{
		((EdgePatchConsSparseSorted*)this)->m_dir = 1 - ((EdgePatchConsSparseSorted*)this)->m_dir;
	}
}

inline TypePartialEnumeration::REAL TypePartialEnumeration::Edge::UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* _buf)
{
	Vector* buf = (Vector*) _buf;
	REAL vMin;

	if (m_type == PatchConsSparseSortedHor)
	{
		int ksource, kdest;
		int mask, i, destLabel, sourceLabel, destInd, sourceInd;
		int* labelsNode1 = Kglobal.m_labels;
		int* sortLabelsNode1 =  Kglobal.m_hor_source;
		int* labelsNode2 = Kglobal.m_labels;
		int* sortLabelsNode2 = Kglobal.m_hor_sink;

		Vector* group_min = (Vector*) _buf;
		REAL* D = (REAL*) _buf;

		int lastpos;

		if (dir == ((EdgePatchConsSparseSorted*)this)->m_dir)
		{
			mask = twopow(Kglobal.m_patchSize*Kglobal.m_patchSize)-1;
			for (i=0;i<Kglobal.m_patchSize;i++)
			{
				mask = mask & (~twopow((i*Kglobal.m_patchSize)));
			}
			for (i=0;i<Kglobal.m_K;i++)
			{
				D[i] = gamma*source->m_data[i] - m_message->m_data[i];
			}
			lastpos = 0;
			i = 0;

			while (i<Kglobal.m_K)
			{
				//Check if source and dest nodes are consistent

				destInd = *(sortLabelsNode2+i);
				destLabel = labelsNode2[destInd];
				sourceInd = *(sortLabelsNode1+lastpos);
				sourceLabel = labelsNode1[sourceInd];

				// Inconsistent labelling
				m_message->m_data[destInd] = std::numeric_limits<TypePartialEnumeration::REAL>::infinity();

				while (((sourceLabel & mask) >> 1) == (destLabel & (mask >> 1))) // find minimum of consistent source labels
				{
					if (D[sourceInd] < m_message->m_data[destInd])
					{
						m_message->m_data[destInd] = D[sourceInd];
					}
					++lastpos;
					if (lastpos == Kglobal.m_K)
					{
						break;
					}
					sourceInd = *(sortLabelsNode1+lastpos);
					sourceLabel = labelsNode1[sourceInd];
				}

				++i;
				if (i < Kglobal.m_K)
				{
					//Fill rest of consistent destination labels
					while  ((destLabel & (mask >>1)) == (labelsNode2[*(sortLabelsNode2+i)]& (mask >>1)))
					{
						m_message->m_data[*(sortLabelsNode2+i)] = m_message->m_data[destInd];
						++i;
						if (i == Kglobal.m_K)
						{
							break;
						}
					}
				}

			}
		}
		else
		{
			mask = twopow(Kglobal.m_patchSize*Kglobal.m_patchSize)-1;
			for (i=0;i<Kglobal.m_patchSize;i++)
			{
				mask = mask & (~twopow((i*Kglobal.m_patchSize)));
			}
			for (i=0;i<Kglobal.m_K;i++)
			{
				D[i] = gamma*source->m_data[i] - m_message->m_data[i];
			}
			lastpos = 0;
			i = 0;
			while (i<Kglobal.m_K)
			{
				//Check if source and dest nodes are consistent
				sourceInd = *(sortLabelsNode2+lastpos);
				sourceLabel = labelsNode2[sourceInd];
				destInd = *(sortLabelsNode1+i);
				destLabel = labelsNode1[destInd];

				// Inconsistent labelling
				m_message->m_data[destInd] = std::numeric_limits<TypePartialEnumeration::REAL>::infinity(); 

				while (((destLabel & mask) >> 1) == (sourceLabel & (mask >> 1))) // find minimum of consistent source labels
				{
					if (D[sourceInd] < m_message->m_data[destInd])
					{
						m_message->m_data[destInd] = D[sourceInd];
					}
					++lastpos;
					if (lastpos == Kglobal.m_K)
					{
						break;
					}
					sourceInd = *(sortLabelsNode2+lastpos);
					sourceLabel = labelsNode2[sourceInd];
				}

				++i;
				if (i < Kglobal.m_K)
				{
					while  (((destLabel & mask) >> 1) == ((labelsNode1[*(sortLabelsNode1+i)]& mask) >> 1))//Fill rest of consistent destination labels
					{
						m_message->m_data[*(sortLabelsNode1+i)] = m_message->m_data[destInd];
						++i;
						if (i == Kglobal.m_K)
						{
							break;
						}
					}
				}
			}
		}

		vMin = m_message->m_data[0];

		for (kdest=1; kdest<Kglobal.m_K; kdest++)
		{
			if (vMin > m_message->m_data[kdest])
							vMin = m_message->m_data[kdest];
		}

		for (kdest=0; kdest<Kglobal.m_K; kdest++)
		{
			m_message->m_data[kdest] -= vMin;
		}
	}
	else if (m_type == PatchConsSparseSortedVert)
	{
		// In the current implementation each super node have the same number of labels.
		// This limitation saves a lot of memory.

		int ksource, kdest;
		int mask,i,destLabel,sourceLabel, destInd, sourceInd;
		int* labelsNode1 = Kglobal.m_labels;
		int* sortLabelsNode1 = Kglobal.m_vert_source;
		int* labelsNode2 = Kglobal.m_labels;
		int* sortLabelsNode2 = Kglobal.m_vert_sink;

		Vector* group_min = (Vector*) _buf;
		REAL* D = (REAL*) _buf;
		REAL intval = 0;
		int lastpos;

		if (dir == ((EdgePatchConsSparseSorted*)this)->m_dir)
		{
			// New Sparse Version
			mask = twopow(Kglobal.m_patchSize*Kglobal.m_patchSize)-1;
			for (i=0;i<Kglobal.m_patchSize;i++)
			{
				mask = mask & (~twopow(i));
			}
			for (i=0;i<Kglobal.m_K;i++)
			{
				D[i] = gamma*source->m_data[i] - m_message->m_data[i];
			}
			lastpos = 0;
			i = 0;

			while (i<Kglobal.m_K)
			{
				//Check if source and dest nodes are consistent
				destInd = *(sortLabelsNode2+i);
				destLabel = labelsNode2[destInd];
				sourceInd = *(sortLabelsNode1+lastpos);
				sourceLabel = labelsNode1[sourceInd];

				// Inconsistent labelling
				m_message->m_data[destInd] = std::numeric_limits<TypePartialEnumeration::REAL>::infinity();

				while (((sourceLabel & mask) >> Kglobal.m_patchSize) == (destLabel & (mask >> Kglobal.m_patchSize))) // find minimum of consistent source labels
				{
					if (D[sourceInd] < m_message->m_data[destInd])
						m_message->m_data[destInd] = D[sourceInd];

					++lastpos;

					if (lastpos == Kglobal.m_K)
						break;

					sourceInd = *(sortLabelsNode1+lastpos);
					sourceLabel = labelsNode1[sourceInd];
				}

				++i;
				if (i < Kglobal.m_K)
				{
					while  ((destLabel & (mask >> Kglobal.m_patchSize)) == (labelsNode2[*(sortLabelsNode2+i)]& (mask >> Kglobal.m_patchSize)))//Fill rest of consistent destination labels
					{
						m_message->m_data[*(sortLabelsNode2+i)] = m_message->m_data[destInd];
						++i;

						if (i == Kglobal.m_K)
							break;
					}
				}
			}		
		}
		else
		{
			mask = twopow(Kglobal.m_patchSize*Kglobal.m_patchSize)-1;
			for (i=0;i<Kglobal.m_patchSize;i++)
			{
				mask = mask & (~twopow(i));
			}

			for (i=0;i<Kglobal.m_K;i++)
			{
				D[i] = gamma*source->m_data[i] - m_message->m_data[i];
			}
			lastpos = 0;
			i = 0;
			while (i<Kglobal.m_K)
			{
				//Check if source and dest nodes are consistent
				sourceInd = *(sortLabelsNode2+lastpos);
				sourceLabel = labelsNode2[sourceInd];
				destInd = *(sortLabelsNode1+i);
				destLabel = labelsNode1[destInd];

				// Inconsistent labelling
				m_message->m_data[destInd] = std::numeric_limits<TypePartialEnumeration::REAL>::infinity(); 

				while (((destLabel & mask) >> Kglobal.m_patchSize) == (sourceLabel & (mask >> Kglobal.m_patchSize))) // find minimum of consistent source labels
				{
					if (D[sourceInd] < m_message->m_data[destInd])
						m_message->m_data[destInd] = D[sourceInd];

					++lastpos;

					if (lastpos == Kglobal.m_K)
						break;

					sourceInd = *(sortLabelsNode2+lastpos);
					sourceLabel = labelsNode2[sourceInd];
				}

				++i;
				if (i < Kglobal.m_K)
				{
					while  (((destLabel & mask) >> Kglobal.m_patchSize) == ((labelsNode1[*(sortLabelsNode1+i)]& mask) >> Kglobal.m_patchSize))//Fill rest of consistent destination labels
					{
						m_message->m_data[*(sortLabelsNode1+i)] = m_message->m_data[destInd];
						++i;
						if (i == Kglobal.m_K)
							break;

					}
				}
			}
		}

		vMin = m_message->m_data[0];

		for (kdest=1; kdest<Kglobal.m_K; kdest++)
		{
			if (vMin > m_message->m_data[kdest])
				vMin = m_message->m_data[kdest];

		}

		for (kdest=0; kdest<Kglobal.m_K; kdest++)
		{
			m_message->m_data[kdest] -= vMin;
		}


	}
	else
	{
		assert(0);
	}

	return vMin;
}

inline void TypePartialEnumeration::Edge::AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir)
{
	assert(ksource>=0 && ksource<Kglobal.m_K);

	int k;

	if (m_type == PatchConsSparseSortedHor)
	{
		int* labelsNode1 = Kglobal.m_labels;
		int* sortLabelsNode1 = Kglobal.m_hor_source;
		int* labelsNode2 = Kglobal.m_labels;
		int* sortLabelsNode2 = Kglobal.m_hor_sink;

		if (dir == ((EdgePatchConsSparseSorted*)this)->m_dir)
		{
			for (k=0; k<Kglobal.m_K; k++)
				dest->m_data[k] += horintcost(labelsNode1[ksource],labelsNode2[k],0, twopow(Kglobal.m_patchSize*Kglobal.m_patchSize),Kglobal.m_patchSize);
		}
		else
		{
			for (k=0; k<Kglobal.m_K; k++)
				dest->m_data[k] += horintcost(labelsNode1[ksource],labelsNode2[k],1, twopow(Kglobal.m_patchSize*Kglobal.m_patchSize),Kglobal.m_patchSize);
		}
	}
	else if (m_type == PatchConsSparseSortedVert)
	{
		int* labelsNode1 =  Kglobal.m_labels;
		int* sortLabelsNode1 = Kglobal.m_vert_source;
		int* labelsNode2 =  Kglobal.m_labels;
		int* sortLabelsNode2 = Kglobal.m_vert_sink;

		if (dir == ((EdgePatchConsSparseSorted*)this)->m_dir)
		{
			for (k=0; k<Kglobal.m_K; k++)
				dest->m_data[k] += vertintcost(labelsNode1[ksource],labelsNode2[k],0, twopow(Kglobal.m_patchSize*Kglobal.m_patchSize),Kglobal.m_patchSize);
		}
		else
		{
			for (k=0; k<Kglobal.m_K; k++)
				dest->m_data[k] += vertintcost(labelsNode1[ksource],labelsNode2[k],1, twopow(Kglobal.m_patchSize*Kglobal.m_patchSize),Kglobal.m_patchSize);
		}
	}

	else
	{
		assert(0);
	}
}

inline TypePartialEnumeration::REAL TypePartialEnumeration::Edge::horintcost(int i, int j, int dir, int K,int patchSize)
{
	int k,mask;

	mask = K-1;
	for (k=0;k<patchSize;k++)
	{
		mask = mask & (~twopow((k*patchSize)));
	}

	if (dir == 0)
	{
		if (((i & mask) >> 1) == (j & (mask >> 1)))
			return 0;
		else
			return std::numeric_limits<TypePartialEnumeration::REAL>::infinity();

	}
	else
	{
		if (((j & mask) >> 1) == (i &( mask >> 1)))
			return 0;
		else
			return std::numeric_limits<TypePartialEnumeration::REAL>::infinity();
	}
}

inline TypePartialEnumeration::REAL TypePartialEnumeration::Edge::vertintcost(int i, int j, int dir, int K, int patchSize)
{
	int k,mask;

	mask = K-1;
	for (k=0;k<patchSize;k++)
	{
		mask = mask & (~twopow(k));
	}

	if (dir == 0)
	{
		if (((i & mask) >> patchSize) == (j & (mask >> patchSize)))
			return 0;
		else
			return std::numeric_limits<TypePartialEnumeration::REAL>::infinity();
	}
	else
	{
		if (((j & mask) >> patchSize) == (i & (mask >> patchSize)))
			return 0;
		else
			return std::numeric_limits<TypePartialEnumeration::REAL>::infinity();
	}
}




//////////////////////////////////////////////////////////////////////////////////

#endif