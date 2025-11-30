#include <stdio.h>

int main()
{
    //sorting algorithm
    int list[] = {2, 4, 8, 1, 2, 5, 8, 3, 2, 3, 5, 7, 7, 4, 2, 1};
    int x, i, j, temp, times, k, a, y, size;


    //Benutze sifeof operator um die Länge des arrays zu bestimmen
    //Der sizeof operator gibt die anzahl an bytes eines datentyps aus bspw:
    //int 2 oder 4 bytes (32 oder 64 bit)
    //float 4 bytes double 8 bytes
    //char 1 byte
    //dividiere dann durch ein element des arrays um die größe nicht in byte sondern in elementen zu erhalten
    
    size = sizeof(list)/sizeof(list[0]); 

    for (i = 0; i<size; i++)
    {
        for (j = 0; j<(size-i); j++)
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

    for (k = 0; k<7; k++)
    {
        a = list[k];
        printf("%d\t", a);
    }
    printf("Die Größe der Liste ist %d", size);
    return 0;
}