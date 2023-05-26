#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

void *print_message_function(void *ptr) {
   char *message;
   message = (char *)ptr;
   printf("%s \n", message);
   return (0);
}

int main() {
   char *message1 = (char *)"Thread 1 is running";
   char *message2 = (char *)"Thread 2 is running";

   pthread_t thread1, thread2;
   int err_2 = pthread_create(&thread2, NULL, print_message_function, (void *)message2);
   int err_1 = pthread_create(&thread1, NULL, print_message_function, (void *)message1);
   pthread_join(thread1, NULL);
   pthread_join(thread2, NULL);

   printf("Thread 1 returns: %d\n", err_1);
   printf("Thread 2 returns: %d\n", err_2);
   exit(0);
}