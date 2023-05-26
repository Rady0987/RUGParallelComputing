#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
int counter = 0;

void *functionC(void *unused)
{
   pthread_mutex_lock(&mutex1);
   counter++;
   printf("Counter value: %d\n", counter);
   pthread_mutex_unlock(&mutex1);
   return NULL;
}

int main()
{
   int err_1, err_2;
   pthread_t thread1, thread2;

   if ((err_1 = pthread_create(&thread1, NULL, &functionC, NULL)))
   {
      printf("Thread creation failed: %d\n", err_1);
   }
   if ((err_2 = pthread_create(&thread2, NULL, &functionC, NULL)))
   {
      printf("Thread creation failed: %d\n", err_2);
   }
   pthread_join(thread1, NULL);
   pthread_join(thread2, NULL);
   printf("Master Thread prints the counter value: %d\n", counter);
   exit(0);
}