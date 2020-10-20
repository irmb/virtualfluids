//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef TBTRANSMITTER_H
#define TBTRANSMITTER_H

#include <string>

/*================================================================================*/
/*  TbTransmitter                                                                 */
/*                                                                                */
/**
This Class provides the base for sending and receiving of data.
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 08.11.07
*/

/*
usage: ...
*/

//////////////////////////////////////////////////////////////////////////
// Transmitter
// macht nichts ausser daten senden und empfangen
template <typename T>
class TbTransmitter
{
public:
    using value_type = T;

public:
    TbTransmitter()          = default;
    virtual ~TbTransmitter() = default;

    virtual bool isLocalTransmitter() const  = 0;
    virtual bool isRemoteTransmitter() const = 0;

    // preprocess (e.g. synchronizing send-/receive-buffer)
    virtual void sendDataSize()    = 0;
    virtual void receiveDataSize() = 0;

    // calculation
    virtual void prepareForSend() {}
    virtual void sendData() = 0;
    virtual void prepareForReceive() {}
    virtual value_type &receiveData() = 0;

    // data-access
    inline value_type &getData() { return this->data; }
    inline const value_type &getData() const { return this->data; }

    // info-section (usable for remote transmitter)
    virtual int getSendToRank() const { return -1; }
    virtual int getSendToTag() const { return -1; }
    virtual int getRecvFromRank() const { return -1; }
    virtual int getRecvFromTag() const { return -1; }

    virtual std::string toString() const = 0;

protected:
    value_type data;
};

#endif // TBTRANSMITTER_H
