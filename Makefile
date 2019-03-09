CFLAGS = -g -Wall -std=c++0x #-Ofast

all: eval eval_fpr eval_load eval_bucket eval_tput

eval: eval.cpp cuckoo.h hashutil.h
	g++ $(CFLAGS) -o eval eval.cpp 

#eval_fpr: eval_fpr.cpp cuckoo.h hashutil.h
#	g++ $(CFLAGS) -o eval_fpr eval_fpr.cpp

#eval_load: eval_load.cpp cuckoo.h hashutil.h
#   g++ $(CFLAGS) -o eval_load eval_load.cpp

#eval_bucket: eval_bucket.cpp cuckoo.h hashutil.h
#	g++ $(CFLAGS) -o eval_bucket eval_bucket.cpp

#eval_tput: eval_tput.cpp cuckoo.h hashutil.h
#	g++ $(CFLAGS) -o eval_tput eval_tput.cpp

clean:
	rm -f eval
	#rm -f eval_fpr
	#rm -f eval_load
	#rm -f eval_bucket
	#rm -f eval_tput
