#include "MWTaskContainer.h"
#include "MWDriver.h"
MWDriver* gimme_the_master();
// Constructor and Destructor
MWTaskContainer::MWTaskContainer()
{
	taskType = MWNORMAL;
	RMC = MWTask::RMC;
	tasks = new MWList<MWTask>("");
}

MWTaskContainer::~MWTaskContainer()
{
	delete tasks;
}

// Task related
void 
MWTaskContainer::pack_work()
{
//	int numtasks = tasks->Number();

//	RMC->pack(&numtasks, 1, 1); //this is done in MWWorker.C
//MWprintf(10,"MWTaskContainer::pack_work\n");
    MWTask * t = tasks->First();
    while ( tasks->AfterEnd() == false ) {
        t = tasks->Current();
		RMC->pack(&(t->number),1,1);
        t->pack_work();
        tasks->Next();
    }
}

void 
MWTaskContainer::unpack_work()
{
//MWDriver *driver = gimme_the_master();
//driver->gimme_a_task();
//MWprintf(10,"MWTaskContainer::unpack_work");
    MWTask * t = tasks->First();
    while ( tasks->AfterEnd() == false ) {
        t = tasks->Current();
		RMC->unpack(&(t->number),1,1);
        t->unpack_work();
        tasks->Next();
    }
}

void
MWTaskContainer::pack_subresults(int tasknum)
{
}

void
MWTaskContainer::unpack_subresults(int tasknum)
{
}

int
MWTaskContainer::FirstNum()
{
	MWTask *t;

	tasks->First();
	t = tasks->Current();
	return t->number;
}

int 
MWTaskContainer::LastNum()
{
	MWTask *t;

    tasks->Last();
    t = tasks->Current();
    return t->number;
}

int
MWTaskContainer::Number()
{
	return tasks->Number();
}

void 
MWTaskContainer::pack_results()
{
		//int results = tasks->Number();

    //RMC->pack(&results, 1, 1);
	MWTask * t = tasks->First();
    while ( tasks->AfterEnd() == false ) {
        t = tasks->Current();
        t->pack_results();
        tasks->Next();
    }
}

void 
MWTaskContainer::unpack_results()
{
		//int results;

    //RMC->unpack(&results, 1, 1);
    MWTask * t = tasks->First();
    while ( tasks->AfterEnd() == false ) {
        t = tasks->Current();
        t->unpack_results();
        tasks->Next();
    }
}

//List related
void 
MWTaskContainer::addTask(MWTask *t)
{
    tasks->Append(t);

}

void
MWTaskContainer::removeAll()
{
	while(tasks->Remove() != NULL){};
}


MWTask* 
MWTaskContainer::First()
{
	return tasks->First();
}

bool 
MWTaskContainer::AfterEnd()
{
	return tasks->AfterEnd();
}

MWTask* 
MWTaskContainer::Current()
{
	return tasks->Current();
}

MWTask* 
MWTaskContainer::Next()
{
	return tasks->Next();
}

void 
MWTaskContainer::printself( int level )
{
    MWprintf ( level, "  Task container from task number %d to %d\n", FirstNum(), LastNum());
}

