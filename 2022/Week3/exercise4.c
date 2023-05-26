#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define LENGTH 5

pthread_mutex_t fill_mutex;
pthread_cond_t cond_var = PTHREAD_COND_INITIALIZER;

int array[LENGTH] = {5, 5, 5, 5, 5};
int flag = 0;

void *fill_array(void *unused)
{
   int i = 0;
   pthread_mutex_lock(&fill_mutex);
   printf("\nEnter values\n");
   for (i = 0; i < LENGTH; i++)
   {
      scanf("%d", &array[i]);
   }
   pthread_mutex_unlock(&fill_mutex);
   pthread_exit(NULL);
}

void *read_values(void *unused)
{
   int i = 0;
   pthread_mutex_lock(&fill_mutex);

   printf("Values filled in array are");
   for (i = 0; i < LENGTH; i++)
   {
      printf(" %4d ", array[i]);
   }
   pthread_mutex_unlock(&fill_mutex);
   pthread_exit(NULL);
}

int main()
{
   pthread_t thread_1, thread_2;
   pthread_attr_t attr;

   int err = pthread_create(&thread_1, NULL, fill_array, NULL);
   err = pthread_create(&thread_2, NULL, read_values, NULL);

   pthread_join(thread_1, NULL);
   pthread_join(thread_2, NULL);

   exit(0);
}
