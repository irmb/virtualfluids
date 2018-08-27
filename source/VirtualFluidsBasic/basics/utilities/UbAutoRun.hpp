#ifndef UB_AUTORUN_HPP
#define UB_AUTORUN_HPP

#define UB_AUTO_RUN(func)                              UB_AUTO_RUN_(func,  __LINE__)
#define UB_AUTO_RUN_(func, nID)                        UB_AUTO_RUN__(func, nID)
#define UB_AUTO_RUN__(func, nID)                       UB_AUTO_RUN___(func, nID)
#define UB_AUTO_RUN___(func, ID)                                                           \
    namespace {                                                                         \
        struct UbAutoRun##ID {                                                            \
            UbAutoRun##ID() {                                                             \
                func;                                                                   \
            }                                                                           \
        } UbAutoRunInst##ID;                                                              \
    }

    // More concise to implement UB_AUTO_RUN using the following, but BCB emits an ICE on it.
    //static bool UB_AutoRun##ID = ( func , false);


#define UB_AUTO_RUN_1(func)                            UB_AUTO_RUN_NAMED(func, 1)                   
#define UB_AUTO_RUN_2(func)                            UB_AUTO_RUN_NAMED(func, 2)                   
#define UB_AUTO_RUN_3(func)                            UB_AUTO_RUN_NAMED(func, 3)                   
#define UB_AUTO_RUN_4(func)                            UB_AUTO_RUN_NAMED(func, 4)                   
#define UB_AUTO_RUN_5(func)                            UB_AUTO_RUN_NAMED(func, 5)                   
                                                       
#define UB_AUTO_RUN_NAMED(func, name)                  UB_AUTO_RUN_NAMED_(func, name, __LINE__)
#define UB_AUTO_RUN_NAMED_(func, name, nID)            UB_AUTO_RUN_NAMED__(func, name, nID)
#define UB_AUTO_RUN_NAMED__(func, name, nID)           UB_AUTO_RUN___(func, _##name##_##nID)
                                                       
#define UB_AUTO_RUN_ONCE(func)                         UB_AUTO_RUN_ONCE_(func,  __LINE__)
#define UB_AUTO_RUN_ONCE_(func, nID)                   UB_AUTO_RUN_ONCE__(func, nID)
#define UB_AUTO_RUN_ONCE__(func, nID)                  UB_AUTO_RUN_ONCE___(func, nID)
#define UB_AUTO_RUN_ONCE___(func, ID)                                                   \
    struct UbAutoRunOnce##ID {                                                            \
        UbAutoRunOnce##ID() {                                                             \
            if (!init()) {                                                              \
                init() = true;                                                          \
                func;                                                                   \
            }                                                                           \
        }                                                                               \
        static bool &init() {                                                           \
            static bool bInit = false;                                                  \
            return bInit;                                                               \
        }                                                                               \
    };                                                                                  \
    static UbAutoRunOnce##ID AutoRunOnceInst##ID;

#define UB_AUTO_RUN_ONCE_1(func)                           UB_AUTO_RUN_ONCE_NAMED(func, 1)                   
#define UB_AUTO_RUN_ONCE_2(func)                           UB_AUTO_RUN_ONCE_NAMED(func, 2)                   
#define UB_AUTO_RUN_ONCE_3(func)                           UB_AUTO_RUN_ONCE_NAMED(func, 3)                   
#define UB_AUTO_RUN_ONCE_4(func)                           UB_AUTO_RUN_ONCE_NAMED(func, 4)                   
#define UB_AUTO_RUN_ONCE_5(func)                           UB_AUTO_RUN_ONCE_NAMED(func, 5)                   
                                                           
#define UB_AUTO_RUN_ONCE_NAMED(func, name)                 UB_AUTO_RUN_ONCE_NAMED_(func, name, __LINE__)
#define UB_AUTO_RUN_ONCE_NAMED_(func, name, nID)           UB_AUTO_RUN_ONCE_NAMED__(func, name, nID)
#define UB_AUTO_RUN_ONCE_NAMED__(func, name, nID)          UB_AUTO_RUN_ONCE___(func, _##name##_##nID)

#endif // ! UB_AUTORUN_HPP
