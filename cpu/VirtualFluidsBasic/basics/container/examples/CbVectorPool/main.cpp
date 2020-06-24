#include "./functions.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
   try
   {
       {
         CbVectorPool<float>* floatPool = new CbVectorPool<float>(0);
         CbVector<float> v1,v2,v3;
         CbVector<float> v0(new CbVectorAllocatorPool<float>(104,floatPool) );
         v0.resize(20);
         v0[3] = 60000;
         v0.resize(40);
         v0[3] = 90000;
         v1.setAllocator( new CbVectorAllocatorPool<float>(100,floatPool) );
         v2.setAllocator( new CbVectorAllocatorPool<float>(101,floatPool) );
         v3.setAllocator( new CbVectorAllocatorPool<float>(102,floatPool) );
         v1.resize(20, 0.0);
         v1.resize(30, 0.0);
         v2.resize(0);
         v2.resize(40, 0.0);
         v3.resize(30, 0.0);
         v3.resize(50, 0.0);

         for(CbVector<float>::size_type i=v1.size()-1; i>=15; i--)
            v1[i] = (CbVector<float>::value_type)i;
         for(CbVector<float>::size_type i=v2.size()-1; i>=35; i--)
            v2[i] = (CbVector<float>::value_type)i;
         for(CbVector<float>::size_type i=v3.size()-1; i>=10; i--)
            v3[i] = (CbVector<float>::value_type)i;
         v1.size(); 
         v2.size();
         v3.size();
         for(CbVector<float>::size_type i=0; i<v1.size(); i++)  v1[i];
         v1.size();
         v2.size();
         v3.size();
         for(CbVector<float>::size_type i=0; i<v2.size(); i++) v2[i];
         v1.size();
         v2.size();
         v3.size();
         for(CbVector<float>::size_type i=0; i<v3.size(); i++) v3[i];
      }
      
     CbVectorPool<value_type>* vectorPool = new CbVectorPool<value_type>(0);

     vector< StlVectorPtr > stlVecs;
     vector< CbVectorPtr >  cbVecs;
     vector< CbVectorPtr >  cbPoolVecs;

     cout<<"check"<<__LINE__<<endl;
     createVecs(10,12,0,stlVecs,cbVecs,cbPoolVecs,vectorPool);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,0,2);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,3,3);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     createVecs(8,5,stlVecs,cbVecs,cbPoolVecs,vectorPool);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,20,7);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,20,3);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,0,7);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,20,3);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     createVecs(4,3,stlVecs,cbVecs,cbPoolVecs,vectorPool);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,20,3);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,0,7);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,20,3);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     //dealloc check
     stlVecs.resize(5);
     cbVecs.resize(5);
     cbPoolVecs.resize(5);

     cout<<"check"<<__LINE__<<endl;
     createVecs(4,3,stlVecs,cbVecs,cbPoolVecs,vectorPool);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     createVecs(4,3,stlVecs,cbVecs,cbPoolVecs,vectorPool);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);


     //operator= check
     CbVector<value_type> testPool1(10, new CbVectorAllocatorPool<value_type>(vectorPool));
     CbVector<value_type> testPool2(1 , new CbVectorAllocatorPool<value_type>(vectorPool));
     CbVector<value_type> testPool3(8 , new CbVectorAllocatorPool<value_type>(vectorPool));
     CbVector<value_type> testPool4(8 , new CbVectorAllocatorPool<value_type>(vectorPool));
     CbVector<value_type> testStd1(10);

     for(CbVector<value_type>::size_type i=0; i<testStd1.size(); i++ )
        testStd1[i] = (value_type)i*10;

     testPool1 = testStd1;
     testPool4 = testStd1;
     testPool3 = testPool4;
     testPool2 = testPool3;

     for(CbVector<value_type>::size_type i=0; i<testStd1.size(); i++ )
        cout<<testStd1[i]<<" "; cout<<endl;
     for(CbVector<value_type>::size_type i=0; i<testPool1.size(); i++ )
        cout<<testPool1[i]<<" "; cout<<endl;
     for(CbVector<value_type>::size_type i=0; i<testPool2.size(); i++ )
        cout<<testPool2[i]<<" "; cout<<endl;
     for(CbVector<value_type>::size_type i=0; i<testPool3.size(); i++ )
        cout<<testPool3[i]<<" "; cout<<endl;
     for(CbVector<value_type>::size_type i=0; i<testPool4.size(); i++ )
        cout<<testPool4[i]<<" "; cout<<endl;
    ///end


     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;
     cout<<"// access test - start"<<endl;
     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;
     cout<<"check"<<__LINE__<<endl;
     createVecs(1000,1000,stlVecs,cbVecs,cbPoolVecs,vectorPool,true);

     CbVectorPool<value_type>* pool2 = new CbVectorPool<value_type>(1);
     vector< StlVectorPtr > stlVecs2;
     vector< CbVectorPtr >  cbVecs2;
     vector< CbVectorPtr >  cbPoolVecs2;
     createVecs(1000,1000,stlVecs2,cbVecs2,cbPoolVecs2,pool2,true);

     cout<<"access check\n";
     //accessCheck(1000,stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     setValues(stlVecs,cbVecs,cbPoolVecs);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"check"<<__LINE__<<endl;
     resize(stlVecs,cbVecs,cbPoolVecs,120,3,true);
     equalCheck(stlVecs,cbVecs,cbPoolVecs);

     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;
     cout<<"// access test - end"<<endl;
     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;

     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;
     cout<<"// EXCEPTION TEST - start"<<endl;
     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;
     delete vectorPool;
     vectorPool = NULL;
     try
     {
        resize(stlVecs,cbVecs,cbPoolVecs,20,3);
     }
     catch(UbException& e)
     {
        cout<<"if exception tells something about \"vectorPool==NULL\" -> test successfully passed:"<<endl;
        cout<<e<<endl;
     }
     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;
     cout<<"// EXCEPTION TEST - end"<<endl;
     cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;

     cout<<"\n\n\nALL TESTS PASSED\n";
   }
   catch(UbException& e)
   {
      std::cerr<<e<<std::endl;
   }
   catch(const std::exception &e)
   {
      std::cerr << "Caught exception:\n";
      std::cerr << "Type: " << typeid(e).name() << "\n";
      std::cerr << "What: " << e.what() << "\n";
   }
   catch(...)
   {
      std::cerr<<" **            Verdammte Scheisse - mal wieder Mist gebaut!                **"<<endl;
   }
    return 0;
}

// #include <functional>
// #include <iostream>
// #include <vector>
// #include <algorithm>
// #include <typeinfo>
//
// struct bar
// {
//    bar()
//       : data(0)
//    {}
//
//    void foo(const std::size_t value) { std::cout << "data = " << value << " (old: " << data << ");" << std::endl; data = value; }
//
// private:
//    std::size_t data;
// };
//
// int main()
// {
//    std::vector<bar> data(10);
//
//    /* operator[] => Indexoperator */
//    for (std::size_t i(0); i < data.size(); ++i)
//       data[i].foo(2);
//
//    /* begin(), end() => Iterator */
//    const std::vector<bar>::iterator it_end(data.end());
//    for (std::vector<bar>::iterator it(data.begin()); it != it_end; ++it)
//       it->foo(3);
//
//    /* for_each => Algorithm | Iterator */
//    std::for_each(data.begin(), data.end(), std::bind2nd(std::mem_fun_ref(&bar::foo), 2));
// }
