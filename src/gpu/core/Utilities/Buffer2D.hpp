#ifndef BUFFER2D_H
#define BUFFER2D_H

#include <iostream>
#include <stdlib.h>

template <class T>
class Buffer2D
{
public:
   Buffer2D()
   {
      data = NULL;
   }
   Buffer2D(int r, int c)
   {
      data = NULL;
      setSize(r, c);
   }
   ~Buffer2D()
   {
      if(data != NULL)
         delete[] data;
   }
   void setSize(int r, int c)
   {
      size = r*c;
      row = r;
      column = c;
      if(data != NULL)
         delete[] data;
      data = new T[size];
   } 
   void Empty()
   {
      if(data != NULL)
      {
         delete[] data;
         data = NULL;
      }
   }
   T& operator [] (const int i)
   {
      try
      {
         if (i > row)
         {
            throw i;
         }
      }
      catch (int i)
      {
         std::cout << "Error: row " << i << " does not exist!" << std::endl;
         exit(EXIT_FAILURE);
      }
      return data[i*column];
   }
   T* getData()
   {
      return data;
   }
   int getSize()
   {
      return size;
   }
   int getRowSize()
   {
      return column;
   }

private:
   T* data;
   int row, column;
   int size;
};
#endif	//BUFFER2D_H
