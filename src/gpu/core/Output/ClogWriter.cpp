#if defined( MPI_LOGGING )

#include "ClogWriter.h"
#include <mpi.h>
#include <mpe.h>

ClogWriter::ClogWriter(int data, char* name, char* color)
{
   //eventID_begin = MPE_Log_get_event_number(); 
   //eventID_end   = MPE_Log_get_event_number(); 
   /*
   user should NOT assign eventIDs directly in MPE_Describe_state()
   Get the eventIDs for user-defined STATES(rectangles) from
   MPE_Log_get_state_eventIDs() instead of the deprecated function
   MPE_Log_get_event_number().
   */
   MPE_Log_get_state_eventIDs(&eventID_begin, &eventID_end);
   MPE_Describe_state(eventID_begin, eventID_end, name, color);
   this->data = data;
}
void ClogWriter::Start()
{
   MPE_Log_event(eventID_begin, data, NULL);
}
void ClogWriter::End()
{
   MPE_Log_event(eventID_end, data, NULL);
}
#endif
