/***************************Copyright-DO-NOT-REMOVE-THIS-LINE**
  *
  * Condor Software Copyright Notice
  * Copyright (C) 1990-2004, Condor Team, Computer Sciences Department,
  * University of Wisconsin-Madison, WI.
  *
  * This source code is covered by the Condor Public License, which can
  * be found in the accompanying LICENSE.TXT file, or online at
  * www.condorproject.org.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * AND THE UNIVERSITY OF WISCONSIN-MADISON "AS IS" AND ANY EXPRESS OR
  * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY, OF SATISFACTORY QUALITY, AND FITNESS
  * FOR A PARTICULAR PURPOSE OR USE ARE DISCLAIMED. THE COPYRIGHT
  * HOLDERS AND CONTRIBUTORS AND THE UNIVERSITY OF WISCONSIN-MADISON
  * MAKE NO MAKE NO REPRESENTATION THAT THE SOFTWARE, MODIFICATIONS,
  * ENHANCEMENTS OR DERIVATIVE WORKS THEREOF, WILL NOT INFRINGE ANY
  * PATENT, COPYRIGHT, TRADEMARK, TRADE SECRET OR OTHER PROPRIETARY
  * RIGHT.
  *
  ****************************Copyright-DO-NOT-REMOVE-THIS-LINE**/
/**************************************************************************
 *  Implementation of MWListElement, MWIndex, MWIndexList, MWList         *
 **************************************************************************/ 

#include "MWList.h"

#if 0

/** Implementation of MWListElement<ObjType> **/

template <class ObjType>
MWListElement<ObjType>::MWListElement(ObjType *itemPtr, double sortKey)
{
	obj= itemPtr;
	key = sortKey;
	next = this; 
	prev = this;
}
	
template <class ObjType>
MWListElement<ObjType>::~MWListElement()
{
	/* Nothing */
}
	
/** Implementation of MWIndex<ObjType> **/
template <class ObjType>
MWIndex<ObjType>::MWIndex<ObjType>()
{
	/* Nothing */
}

template <class ObjType>
MWIndex<ObjType>::~MWIndex<ObjType>()
{
	/* Nothing */
}

template <class ObjType>
MWIndex<ObjType>*
MWIndex<ObjType>::Split()
{
	int i;

	if (this->num<2)
		return NULL;
	
	MWListElement<ObjType> * current = this->first->prev;

	/* Walk through the list until we meet the half point */
	for (i=0; i<num/2; i++)
		current = current->next;
	
	/* Create a new index for the second half */
	MWIndex<ObjType> * i2 = new MWIndex<ObjType>;
	if (i2 == NULL)
		return NULL;

	i2->num = num - i;
	i2->first = current->next;
	i2->lower = current->key; /* must not leave space between boundaries */
	i2->last = this->last;
	i2->upper = this->upper;
	
	/* Update the current index */
	this->num = num/2;
	this->last = current;
	this->upper = current->key;

	/* return the new index, which should be inserted into index list */
	return i2;
}

/* split an index into many smaller ones whose sizes <= bucket_size */
template <class ObjType>
void
MWIndex<ObjType>::SplitTo(MWIndexList<ObjType>* ilist, int bucket_size)
{
	if (num < bucket_size)
		return;
	else {
		ilist->Split(this);
		this->SplitTo(ilist, bucket_size);
		this->next->SplitTo(ilist, bucket_size);
	}
}

/* Balance the current index with the next index, either split or merge, 
 * to make small buckets together, and cut large buckets into halves. 
 * Return ptr to new index (current when splitting or next when merging) */
template <class ObjType>
MWIndex<ObjType>*
MWIndex<ObjType>::BalanceNext(MWIndexList<ObjType>* ilist, int bucket_size)
{
	/* Need to split */
	if (this->num > bucket_size) {
		MWIndex<ObjType> * next = this->next;
		this->SplitTo(ilist, bucket_size+1);
		return next;
	}

	/* Already the end of index list */
	if (this->next == ilist->dummy)
		return this->next;

	/* Need to merget */
	MWIndex<ObjType>* next_ind = this->next;
	if (this->num + next_ind->num < bucket_size+1) {
		return ilist->Merge(this);
	} else {  
		/* Nothing to be changed */
		return this->next;
	}
}

/** Implementation of MWIndexList<ObjType> **/
template <class ObjType>
MWIndexList<ObjType>::MWIndexList()
{
	dummy = new MWIndex<ObjType>();
	if (dummy == NULL) {
		MWprintf(10, "MWList::Can't initialize the list during object construction.\n");
		return;
	}
	dummy->num = 0;
	dummy->next = dummy;
	dummy->prev = dummy;
	dummy->lower = -DBL_MAX;
	dummy->upper = DBL_MAX;
	dummy->first = NULL;
	dummy->last = NULL;

	num = 0;
}

template <class ObjType>
MWIndexList<ObjType>::~MWIndexList()
{
	MWIndex<ObjType>* ptr = dummy;
	
	while (num !=0) {
		ptr = ptr->next;
		delete ptr;
		num --;
	}
	delete dummy;

	MWprintf(10, "MWIndexList destructed.\n");
}

/* Build index for the list, split index whose number >= bucket_size,
 * and merge adjacent index whose number < 0.25 * bucket_size */
template <class ObjType>
void 
MWIndexList<ObjType>::BuildIndex(MWList<ObjType>* list, int bucket_size)
{ 
	MWprintf(10, "MWIndexList::BuildIndex()\n");

	if (list->Number() < BUCKET_SIZE/2) {
		MWprintf(10, "MWIndexList:: \
			Will not build index when the list is still short.\n");
		return;
	}
	
	if (num == 0) { /* create index */
		list->Sort();
		
		/* create the first index for the whole list */
		MWIndex<ObjType> * ind = new MWIndex<ObjType>;
		ind->num = list->Number();
		ind->upper = DBL_MAX;
		ind->lower = -DBL_MAX;
		ind->first = list->dummy->next;
		ind->last= list->dummy->prev;

		/* Insert into index list */
		dummy->next = ind;
		ind->next = dummy;
		dummy->prev = ind;
		ind->prev = dummy;

		num ++;
		
		/* split into small enough buckets */
		ind->SplitTo(this, bucket_size);
		
	} else { 
		/* balance index items */
		MWIndex<ObjType> * ind = dummy->next;
		
		while ( (ind != this->dummy) && (ind != NULL) ) {
			ind = ind->BalanceNext(this, bucket_size);
		}
	}
}

template <class ObjType>
void
MWIndexList<ObjType>::Split(MWIndex<ObjType>* index)
{ 
	MWIndex<ObjType> * another = index->Split();
	index->next->prev = another;
	another->next = index->next;
	index->next = another;
	another->prev = index;

	num ++;
}
	
template <class ObjType>
MWIndex<ObjType>* 
MWIndexList<ObjType>::IndexFor(double key)
{
	MWIndex<ObjType> * ptr = dummy;
	for (int i=0; i<num; i++) {
		ptr = ptr->next;
		if ( (ptr->upper >= key) && (ptr->lower <= key) )
			return ptr;
	}
	return NULL;
}

/* Merge two index, and delete any empty index */
template <class ObjType>
MWIndex<ObjType>*
MWIndexList<ObjType>::Merge(MWIndex<ObjType> * i1)
{
	if (num == 1) {/* Can't merge */
		if (i1->num == 0) { /* Delete i1 then */
			dummy->next = dummy->prev = dummy;
			num = 0;
			delete i1;
			return NULL;
		}
	}
	
	if (i1->next != dummy) {
		/* need to merge with it's next */
	       	MWIndex<ObjType> * next_ind = i1->next;
	       	i1->num = i1->num + next_ind->num;
	       	i1->last = next_ind->last;
	       	i1->upper = next_ind->upper;

		/* remove the item from the list */
	       	next_ind->next->prev = i1;
	       	i1->next = next_ind->next;
	       	delete next_ind;
	       	num --;
	}

	return i1;
}

template <class ObjType>
void
MWIndexList<ObjType>::Display()
{
	MWIndex<ObjType> * ind = dummy;
	int total = 0;
	for (int i=0; i<num; i++) {
		ind = ind->next;
		if (ind == NULL) {
			MWprintf(10, "MWList:: Bad problem, index list broken.\n");
		}
			
		total += ind->num;
		if (ind->num > 0) {
			MWprintf(10, "[%d] num=%4d lower=%12.2f\tupper=%12.2f\tfirst=%12.2f\tlast=%12.2f\n", 
				i, ind->num, ind->lower, ind->upper, ind->first->key, ind->last->key);
		}
	}
	MWprintf(10, "Total num = %d.\n", total);
}

/** Implementation of MWList<ObjType> **/

/* The empty list contains a dummy element which has a null ptr 
 * to an object. dummy->next will be the first obj, dummy->prev 
 * will be the last obj when the list is not empty. */
template <class ObjType>
MWList<ObjType>::MWList()
{
	dummy = new MWListElement<ObjType>(NULL, -DBL_MAX);
	if (dummy == NULL) {
		MWprintf(10, "MWList::Can't initialize the list during object construction.\n");
		return;
	}
	dummy->next = dummy;
	dummy->prev = dummy;
	
	current = dummy;
	num = 0;

	/* Indexing */
	sorted = false;
	index = NULL;
}

template <class ObjType>
MWList<ObjType>::~MWList()
{ 
	if (num == 0)
		return;
	current = dummy->next;
	while (RemoveCurrent() != NULL); // delete all the list elements
	delete dummy;
}

template <class ObjType>
bool
MWList<ObjType>::IsEmpty() 
{
       return (num == 0);
}

template <class ObjType>
int
MWList<ObjType>::Number()
{
	return num;
}

template <class ObjType>
ObjType*
MWList<ObjType>::First()
{
	if (num == 0)
		return NULL;
	
	current = dummy->next;
	return current->obj;
}

template <class ObjType>
ObjType*
MWList<ObjType>::Last()
{
	if (num == 0)
		return NULL;
	
	current = dummy->prev;
	return current->obj;
}

template <class ObjType>
bool 
MWList<ObjType>::BeforeHead()
{
	return (current == dummy);
}

template <class ObjType>
bool 
MWList<ObjType>::AfterEnd()
{
	return (current == dummy);
}

template <class ObjType>
void
MWList<ObjType>::ToHead()
{
	current = dummy->next;
}

template <class ObjType>
void
MWList<ObjType>::ToEnd()
{
	current = dummy->prev;
}

/* Return the current object, unless currently at head or end */
template <class ObjType>
ObjType *
MWList<ObjType>::Current()
{
	if (num==0)
		return NULL;
	
	if (current == dummy)
		return NULL;

	return current->obj;
}

/* Return the position, can be used to remove a particular element */
template <class ObjType>
MWListElement<ObjType> *
MWList<ObjType>::CurPosition()
{
	return current;
}

/* Advance the current ptr, and return the new obj that current pointed to */
template <class ObjType>
ObjType *
MWList<ObjType>::Next()
{
	current = current->next;
	if (current== dummy)
		return NULL;
	else return current->obj;
}

/* Retreat the current ptr, and return the new obj that current pointed to */
template <class ObjType>
ObjType *
MWList<ObjType>::Prev()
{
	current = current->prev;
	if (current->prev == dummy)
		return NULL;
	else return current->obj;
}

/* Remove an element from the list */
template <class ObjType>
void 
MWList<ObjType>::Remove(MWListElement<ObjType> * item)
{
	if ((item == dummy)||(num == 0)) {
		MWprintf(10, "MWList::Can't remove the current element from empty list.\n");
		return;
	}
	
	/* Update index */
	if (index != NULL) {
	       	MWIndex<ObjType> * ind = index->IndexFor(item->key);
	       	if (ind->num == 1)
		       	index->Merge(ind);

		if (ind->num > 1) {
		       	if (item == ind->last) 
				ind->last = item->prev; 
			if (item == ind->first)
			       	ind->first = item->next;
	       	} else {
		       	ind->first = ind->last = NULL;
	       	}

		ind->num --;
	}
	
	/* Remove the item */
	item->prev->next = item->next;
	item->next->prev = item->prev;
	delete item;
	num --;
	
	if ((num == 0))
		ClearIndex();
}

/* Remove the current element, and return the obj it points to. 
 * The current pointer will become the previous element, or if removing
 * the first element, the dummy position */
template <class ObjType>
ObjType *
MWList<ObjType>::RemoveCurrent()
{
	if (num == 0) {
		MWprintf(10, "MWList::Can't remove any element from empty list.\n");
		return NULL;
	}
	
	ObjType * newObj = current->obj;
	current = current->prev;
	Remove(current->next);
	return newObj;
}

/* Rewind to head position, remove the first element, 
 * and return the object it points to. - to keep compatible with the old API*/
template <class ObjType>
ObjType *
MWList<ObjType>::Remove()
{
	if (num == 0) {
		MWprintf(10, "MWList::Can't remove any element from empty list.\n");
		return NULL;
	}
	
	ObjType * newObj = dummy->next->obj;
	Remove(dummy->next);
	current = dummy;
	return newObj;
}

/* Pop the first element from the list */
template <class ObjType>
ObjType*
MWList<ObjType>::Pop()
{
	return Remove();
}

/* Append an object at the end of the list */
template <class ObjType>
MWListElement<ObjType>*
MWList<ObjType>::Append(ObjType *obj)
{
	MWListElement<ObjType> * item = new MWListElement<ObjType>(obj, 0.0);
	
	if (item == NULL) 
		return NULL;
	
	dummy->prev->next = item;
	item->prev = dummy->prev;
	dummy->prev = item;
	item->next = dummy;
	current = item;
	num ++;

	if (num == 1)
		sorted = true;

	if (sorted)  { /* Update sorted */
		if ( (num>1) && (dummy->prev->prev->key > 0.0) ) {
			sorted = false;
			if (index != NULL) 
				ClearIndex();
		} else if (index != NULL) { /* Update index */
		       	MWIndex<ObjType> * ind = index->IndexFor(0.0);
		       	if (ind->num == 1)
			       	index->Merge(ind);

			if (ind->num > 1) {
			       	if (item == ind->last) 
					ind->last = item->prev; 
				if (item == ind->first)
				       	ind->first = item->next;
		       	} else {
			       	ind->first = ind->last = NULL;
		       	}

			ind->num --;
		}
	}
	
	return item;
}

/* Prepend an object at the head of the list */
template <class ObjType>
MWListElement<ObjType>*
MWList<ObjType>::Prepend(ObjType *obj)
{
	MWListElement<ObjType> * item = new MWListElement<ObjType>(obj, 0.0);
	
	if (item == NULL) 
		return NULL;
	
	dummy->next->prev = item;
	item->next = dummy->next;
	dummy->next = item;
	item->prev = dummy;
	current = item;
	num ++;

	if (num == 1)
		sorted = true;

	if (sorted) { /* Update sorted */
		if ( (num>1) && (dummy->next->next->key < 0.0) ) {
			sorted = false;
			MWprintf(10, "MWList becomes not sorted.\n");
			if (index != NULL)
				ClearIndex();
		} else if (index != NULL) { /* Update index */
		       	MWIndex<ObjType> * ind = index->IndexFor(0.0);
		       	if (ind->num == 1)
			       	index->Merge(ind);

			if (ind->num > 1) {
			       	if (item == ind->last) 
					ind->last = item->prev; 
				if (item == ind->first)
				       	ind->first = item->next;
		       	} else {
			       	ind->first = ind->last = NULL;
		       	}

			ind->num --;
		}
	}

	return item;
}

template <class ObjType>
void 
MWList<ObjType>::Display()
{
	MWprintf(80, "Number = %d\n", num);
#if 0
	MWListElement<ObjType> *p;
	
	MWprintf(80, "<< -- \n");
       	p = dummy->next;
	while (p!=dummy){
		MWprintf(80, "  %f\t", p->key);
		p->obj->Display();
		p = p->next;
	}
	MWprintf(80, "-- >>\n");
#endif
}

/* Add an object into the list in a sorted order */
template <class ObjType>
void 
MWList<ObjType>::SortedInsert(ObjType *obj, double key)
{
	MWListElement<ObjType> * ptr;

	if ( (sorted == false) && (num > 0) ) {
		MWprintf(10, "MWList::SortedInsert into unsorted list.\n");
		exit(1);
	}
		
	if (index != NULL) { /* indexed */
		IndexedInsert(obj, key);
		return;
	} 

	/* If not indexed */
	if (num == 0) {
		MWListElement<ObjType>* ele = new MWListElement<ObjType>(obj, key);
		dummy->next = dummy->prev = ele;
		ele->next = ele->prev = dummy;
		num ++;
		sorted = true;
	} else { 
		/* walk through the list to find a place to insert */
		ptr = dummy->next;
		while ( (ptr != dummy) && (key > ptr->key) ) 
				ptr = ptr->next;
		
		MWListElement<ObjType> * item = 
			new MWListElement<ObjType>(obj, key);
		/* Insert before ptr */
		ptr->prev->next = item;
		item->prev = ptr->prev;
		item->next = ptr; 
		ptr->prev = item;
		num ++;

		/* Build index when there are already many elements */
		if ( (num > BUCKET_SIZE) && (index == NULL) )
			this->BuildIndex();
	}
}

/* Delete the first item with given key from the sorted list */
template <class ObjType>
ObjType *
MWList<ObjType>::SortedRemove(double key)
{
	MWListElement<ObjType> * ptr;

	if ( (sorted == false) && (num > 0) ) {
		MWprintf(10, "MWList::SortedRemove from unsorted list.\n");
		return NULL;
	}
		
	if (index != NULL) /* indexed */
		return IndexedRemove(key);

	/* If not indexed */
	if (num == 0) {
		return NULL;
	} else { 
		/* walk through the list to find the element with given key */
		ptr = dummy->next;
		while ( (ptr != dummy) && (key > ptr->key) ) 
				ptr = ptr->next;
		
		if ( (key < ptr->key) || (ptr == dummy) ) /* not found */
			return NULL;
		
		/* delete ptr*/
		ObjType * obj = ptr->obj;
		
		ptr->prev->next = ptr->next;
		ptr->next->prev = ptr->prev;
		delete ptr;
		num --;

		return obj;
	}
}

/* With the help of index, add an object into an ordered list */
template <class ObjType>
void 
MWList<ObjType>::IndexedInsert(ObjType *obj, double key)
{  
	static int newly_inserted = 0;

	MWIndex<ObjType> * ind = index->IndexFor(key);
	if (ind == NULL) {
		MWprintf(10, "MWList::IndexedInsert() can't find an index for key = %f.\n", key);
		return;
	}
       	/* walk through the list to find a place to insert */
       	MWListElement<ObjType>* ptr = ind->first;
       	while ( (ptr != dummy) && (key > ptr->key) ) 
		ptr = ptr->next;

	MWListElement<ObjType> * item = 
		new MWListElement<ObjType>(obj, key);
	
       	/* Insert before ptr */
       	ptr->prev->next = item;
       	item->prev = ptr->prev;
       	item->next = ptr; 
	ptr->prev = item;
	ind->num ++;
       	num ++;
	
	/* Need to update the index */
	if (ptr == ind->last->next) 
		ind->last = item;
	if (ptr == ind->first) 
		ind->first = item;

	/* Balance */
	newly_inserted ++;
	
	if (newly_inserted > BUCKET_SIZE) {
		BuildIndex();
		newly_inserted = 0;
	}

}

/* With the help of index, add an object into an ordered list */
template <class ObjType>
ObjType *
MWList<ObjType>::IndexedRemove(double key)
{  
	static int newly_removed = 0;
	
	if (num == 0) 
		return NULL;
	
	MWIndex<ObjType> * ind = index->IndexFor(key);
	if (ind == NULL) {
		MWprintf(10, "MWList::IndexedRemove() can't find an index.\n");
		return NULL;
	}
	
	/* walk through the list to find the element with given key */
       	MWListElement<ObjType> * ptr = ind->first;
       	while ( (ptr != dummy) && (key > ptr->key) ) 
		ptr = ptr->next;

	if ( (key < ptr->key) || (ptr == dummy) ) /* not found */
	       	return NULL;

	/* Need to update the index */
	if (ind->num == 1)
		index->Merge(ind);
	
	if (ind->num > 1) {
		if (ptr == ind->last) 
			ind->last = ptr->prev; 
		if (ptr == ind->first)
		       	ind->first = ptr->next;
	} else {
		ind->first = ind->last = NULL;
	}
	
	/* delete ptr*/
       	ObjType * obj = ptr->obj;
	ptr->prev->next = ptr->next;
       	ptr->next->prev = ptr->prev;
       	delete ptr;
	ind->num --;
       	num --;

	/* Balance */
	if (num == 0) 
		ClearIndex(); 
	else {
		newly_removed ++;
		if (newly_removed > BUCKET_SIZE) {
			BuildIndex();
			newly_removed = 0;
		}
	}
	
	return obj;
}

/* Sort the list */
template <class ObjType>
void 
MWList<ObjType>::Sort()
{
	if (sorted == true)
		return;

	MWList<ObjType> * list = new MWList<ObjType>;
	ObjType * obj;
	double key;
	
	int inserted = 0;
	while (num > 0) {
		obj = dummy->next->obj;
		key = dummy->next->key;
		list->SortedInsert(obj, key);
		Remove(dummy->next);
		inserted ++;

		/* Use index to accelerate insertion */
		if ( (index == NULL)&&(inserted > BUCKET_SIZE/2) ) {
			list->BuildIndex();
			inserted = 0;
		}
	}

	sorted = true;
}

/* Build index */
template <class ObjType>
void 
MWList<ObjType>::BuildIndex()
{  
	if (index == NULL) 
		index = new MWIndexList<ObjType>;

	index->BuildIndex(this);
	index->Display();
}

/* Build index */
template <class ObjType>
void 
MWList<ObjType>::ClearIndex()
{  
	if (index == NULL) 
		return;

	delete index;
	index = NULL;
}

#endif
