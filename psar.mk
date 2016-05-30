OBJS = psar.o model.o sequence.o scores.o params.o
EXEC = psar
CXX = g++ 

$(EXEC) : $(OBJS)
	$(CXX) -o $@ $^

$(OBJS) : sequence.h
psar.o : model.h scores.h params.h
model.o : model.h params.h
scores.o : scores.h 
params.o : params.h 

.PHONY : clean
clean :
	rm $(EXEC) $(OBJS)
