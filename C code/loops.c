#include <stdio.h>

int main()
{
    //sorting algorithm
    int x, i, j, temp, times;
    int list[] = {4, 1, 2, 8, 3, 14, 10, 9};
    for (i = 0; i<7; i++)
    {
        for (j = 0; j<(7-i); j++)
        {
        x = list[j] - list[j+1]; //Wenn das i-te Element größer als das i+1-te Element ist müssen die plätze getauscht werden
        if (list[j] > list[j+1])
        {
            temp = list[j+1];
            list[j+1] = list[j];
            list[j] = temp;

        }
        }
        times++;
    }

    int k, a;
    for (k = 0; k<7; k++)
    {
        a = list[k];
        printf("%d\t", a);
    }
    return 0;
}