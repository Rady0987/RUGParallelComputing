#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define LENGTH 1000
#define THREAD_COUNT 100

int vector[LENGTH];

void merge(int low, int mid, int high)
{
    int *left = malloc(sizeof(*left) * (mid - low + 1));
    int *right = malloc(sizeof(*right) * (high - mid));
    int length_left = mid - low + 1, length_right = high - mid;

    for (int i = 0; i < length_left; i++)
        left[i] = vector[i + low];
    for (int i = 0; i < length_right; i++)
        right[i] = vector[i + mid + 1];

    int i = 0, j = 0, k = low;
    while (i < length_left && j < length_right)
    {
        if (left[i] <= right[j])
            vector[k++] = left[i++];
        else
            vector[k++] = right[j++];
    }

    while (i < length_left)
    {
        vector[k++] = left[i++];
    }

    while (j < length_right)
    {
        vector[k++] = right[j++];
    }
}

void merge_sort(int low, int high)
{
    int mid = (low + high) / 2;
    if (low < high)
    {
        merge_sort(low, mid);
        merge_sort(mid + 1, high);
        merge(low, mid, high);
    }
}

void *merge_sort_thread(void *arg)
{
    unsigned long thread_part = (unsigned long)arg;
    int remainder = LENGTH % THREAD_COUNT;
    int low = thread_part * (LENGTH / THREAD_COUNT);
    int high = (thread_part + 1) * (LENGTH / THREAD_COUNT) - 1;

    if (thread_part == THREAD_COUNT - 1)
        high += remainder;

    int mid = low + (high - low) / 2;

    if (low < high)
    {
        merge_sort(low, mid);
        merge_sort(mid + 1, high);
        merge(low, mid, high);
    }
    return 0;
}

void merge_all_sections(int numberThreads, int depthFactor)
{
    int valPerThread = LENGTH / THREAD_COUNT;

    for (int thread_part = 0; thread_part < numberThreads; thread_part += 2)
    {
        int low = thread_part * valPerThread * depthFactor;
        int high = ((thread_part + 2) * valPerThread * depthFactor) - 1;
        int mid = low + (valPerThread * depthFactor) - 1;

        if (high >= LENGTH)
            high = LENGTH - 1;

        merge(low, mid, high);
    }

    if (numberThreads / 2 >= 1)
        merge_all_sections(numberThreads / 2, depthFactor * 2);
}

int main()
{
    pthread_t numberThreads[THREAD_COUNT];
    for (int i = 0; i < LENGTH; i++)
        vector[i] = rand() % 100;
    for (int i = 0; i < THREAD_COUNT; i++)
        pthread_create(&numberThreads[i], NULL, merge_sort_thread, (void *)i);
    for (int i = 0; i < THREAD_COUNT; i++)
        pthread_join(numberThreads[i], NULL);

    merge_all_sections(THREAD_COUNT, 1);

    for(int i = 0; i< LENGTH; i++) {
        printf("%d ", vector[i]);
    }

    return 0;
}