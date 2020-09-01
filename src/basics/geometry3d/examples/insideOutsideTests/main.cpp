#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbSphere3D.h>
#include <geometry3d/GbCylinder3D.h>

int main(int argc, char** argv)
{
     GbSphere3D test(10,10,10,8);
     
     bool cutSp1    /*false*/= test.isCellCuttingGbObject3D(9,9,9,11,11,11); //cell komplett IN sphere
     bool cutSp2    /*true */= test.isCellCuttingGbObject3D(0,0,0,20,20,20); //cell umhuellt    sphere
     bool cutSp3    /*true */= test.isCellCuttingGbObject3D(0,0,0,10,10,10); //cell cutted      sphere
     bool cutSp4    /*false*/= test.isCellCuttingGbObject3D(100,100,100,101,101,101); //cell nix      sphere

     
     bool cutInsSp1 /*true */= test.isCellInsideOrCuttingGbObject3D(9,9,9,11,11,11);     //cell komplett IN sphere
     bool cutInsSp2 /*true */= test.isCellInsideOrCuttingGbObject3D(0,0,0,20,20,20);     //cell umhuellt    sphere
     bool cutInsSp3 /*true */= test.isCellInsideOrCuttingGbObject3D(0,0,0,10,10,10);     //cell cutted      sphere
     bool cutInsSp4 /*false*/= test.isCellInsideOrCuttingGbObject3D(100,100,100,101,101,101); //cell nix      sphere
                                                                    
     GbCuboid3D test1(0,0,0,10,10,10);

     bool cutCu1    /*false*/= test1.isCellCuttingGbObject3D(4,4,4,6,6,6);       //cell komplett IN cube
     bool cutCu2    /*true */= test1.isCellCuttingGbObject3D(-1,-1,-1,11,11,11); //cell umhuellt    cube
     bool cutCu3    /*true */= test1.isCellCuttingGbObject3D(5,5,5,15,15,15);   //cell cutted      cube
     bool cutCu4    /*false*/= test1.isCellCuttingGbObject3D(12,12,12,15,15,15);   //cell nix      cube

     bool cutInsCu1 /*true */= test1.isCellInsideOrCuttingGbObject3D(4,4,4,6,6,6);       //cell komplett IN cube
     bool cutInsCu2 /*true */= test1.isCellInsideOrCuttingGbObject3D(-1,-1,-1,11,11,11); //cell umhuellt    cube
     bool cutInsCu3 /*true */= test1.isCellInsideOrCuttingGbObject3D(5,5,5,15,15,15);   //cell cutted      cube
     bool cutInsCu4 /*false*/= test1.isCellInsideOrCuttingGbObject3D(12,12,12,15,15,15);  //cell nix      cube

     GbCylinder3D test2( 0,0,0, 20, 0, 0, 10);

     bool cutCy1     /*false*/ = test2.isCellCuttingGbObject3D(1,-1,-1,4,1,1);       //cell komplett IN cyl
     bool cutCy2     /*true */ = test2.isCellCuttingGbObject3D(10,0,0,15,12,11);     //cell umhuellt    cyl
     bool cutCy3a    /*true */ = test2.isCellCuttingGbObject3D(5,5,5,15,15,15);      //cell cutted      cyl im kreisbreich
     bool cutCy3b    /*true */ = test2.isCellCuttingGbObject3D(-5,-1,-1,5,1,1);      //cell cutted      cyl am stirn
     bool cutCy4     /*false*/= test2.isCellCuttingGbObject3D(-10,-10,-10,-5,-5,-5);   //cell nix      cyl


     bool cutInsCy1  /*true */= test2.isCellInsideOrCuttingGbObject3D(4,4,4,6,6,6);       //cell komplett IN cube
     bool cutInsCy2  /*true */= test2.isCellInsideOrCuttingGbObject3D(10,0,0,15,12,11);   //cell umhuellt    cyl
     bool cutInsCy3a /*true */= test2.isCellInsideOrCuttingGbObject3D(5,5,5,15,15,15);      //cell cutted      cyl im kreisbreich
     bool cutInsCy3b /*true */= test2.isCellInsideOrCuttingGbObject3D(-5,-1,-1,5,1,1);   //cell cutted      cube
     bool cutInsCy4  /*false*/= test2.isCellInsideOrCuttingGbObject3D(-10,-10,-10,-5,-5,-5);   //cell nix      cube
}