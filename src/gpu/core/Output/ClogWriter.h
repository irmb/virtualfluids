#ifndef CLOGWRITER_H
#define CLOGWRITER_H

#if defined( MPI_LOGGING )

class ClogWriter
{
public:
   ClogWriter(int data, char* name, char* color);
   void Start();
   void End();
protected:
private:
   int eventID_begin, eventID_end;
   int data;
};
#endif
#endif    //CLOGWRITER_H
