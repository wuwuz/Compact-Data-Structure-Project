CFLAGS = -g -Wall -std=c++0x

all: eval eval_item eval_load

eval: eval.cpp cuckoo.h hashutil.h
	g++ $(CFLAGS) -o eval eval.cpp 

eval_item: eval_item.cpp cuckoo.h hashutil.h
	g++ $(CFLAGS) -o eval_item eval_item.cpp

eval_load: eval_load.cpp cuckoo.h hashutil.h
	g++ $(CFLAGS) -o eval_load eval_load.cpp

clean:
	rm -f eval
	rm -f eval_item
	rm -r eval_load
