/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef NON_CREATABLE_H
#define NON_CREATABLE_H

class NonCreatable
{
private:
    NonCreatable()                     = delete;
    NonCreatable(const NonCreatable &) = delete;
    NonCreatable &operator=(const NonCreatable &) = delete;
};

#endif
