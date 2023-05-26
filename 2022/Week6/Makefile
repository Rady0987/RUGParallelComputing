# Set the compilation settings
CC=mpicc
CFLAGS=-Wall -Wextra -pthread -O2
LIBS=-lcrypt -lpthread

# Debug settings
ifdef DEBUG
CFLAGS += -g
endif

# Set where to leave our compiled intermediate object files
ODIR=obj

# List of libraries
_OBJS = hashtable.o linkedlist.o pw_helpers.o stringlist.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS)) # Make prepended with obj/

# Name of our implementation
IMPL = guessword

#########################

all: $(IMPL)

# Rule to automatically compile any object file
$(ODIR)/%.o: %.c %.h
	@mkdir -p $(@D)
	@$(CC) -c -o $@ $< $(CFLAGS)

# Rule to build reference implementation
$(IMPL): %: $(OBJS) %.c
	@$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# Cleanup
.PHONY: clean # Prevent make from confusing this with a file
clean:
	@rm -rf $(ODIR)/*.o $(ODIR)
	@rm -f $(IMPL)
	$(info Cleaned up)
