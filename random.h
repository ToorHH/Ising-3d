#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED
#include <math.h>
const int rand_multiplier = 16807;
const int rand_modulus = 0x7fffffff;//2^31-1
const int rand_quotient = 127773;   //q=m/a;
const int rand_remainder = 2836;    //r=m%a;
const int rand_range = 0x7ffffffe;  //MODULOUS-1

//Last Modified 2015/01/06
/*
*16807随机数生成器,生成[0,2^31-1]内均匀分布的随机数
*种子值为CPU时钟周期数TSC（Time Stamp Counter）
*用线性同余法产生下一个值new_seed=(a*seed+b) mod m
*/
class rand_16807
{
public:
    rand_16807()                    //初始种子值为TSC除以2
    :   __init_seed(__get_cycle_count() >> 1),
        __seed(__init_seed)
    {}

    rand_16807(int init_seed)       //用指定种子值初始化
    :   __init_seed(init_seed),
        __seed(__init_seed)
    {}
    int seed(int new_seed = -1)     //实现了种子的变更和重置
    {
        if (new_seed == -1)
            __seed = __init_seed;                 //当没有传入新值，将种子值重置为初始值
        else
            __seed = new_seed;                    //当传入了新的种子值便更新
        return __seed;                            //返回当前种子值
    }
    int operator()()
    {
        __seed = rand_multiplier * (__seed % rand_quotient) - rand_remainder * (__seed / rand_quotient);//schrage法取模
        if (__seed < 0)
            __seed += rand_modulus;
        return __seed;              //前进一次并返回产生的随机数
    }
    int advance(int n)              //前进n步并返回产生的随机数
    {
        while (n-- > 0)
            operator()();
        return __seed;
    }
private:
    unsigned int __get_cycle_count()
    {
        unsigned int ret;
        __asm__ ("RDTSC":"=a"(ret));
        return ret;
    }
    int __init_seed;
    int __seed;
};
class uni_real_dis{     //由整数的随机数产生器rand_16807产生均匀分布
public:
    uni_real_dis()      //默认取值区间为[0,1]
    :   __min_x(0),
        __max_x(1),
        __range(1)
    {}
    uni_real_dis(double min_x, double max_x)    //初始化时将取值区间设为[min_x,max_x]
    :   __min_x(min_x),
        __max_x(max_x),
        __range(max_x - min_x)
    {}
    double operator()(rand_16807& rand)         //将rand_16807的输出乘以均匀分布的区间长度再除以2^31-1最后加上偏移量min_x
    {
        return rand() * __range / rand_range + __min_x;
    }
    double set_range(double min_x, double max_x)//将取值区间设为[min_x,max_x]
    {
        __min_x = min_x;
        __max_x = max_x;
        return (__range = max_x - min_x);
    }
private:
    double __min_x;
    double __max_x;
    double __range;
};
class gauss_real_dis{           //由整数的随机数产生器rand_16807产生正态分布
public:
    gauss_real_dis()            //默认参数为0,1
    :   __mu(0),
        __sigma(1)
    {}
    gauss_real_dis(double mu, double sigma) //初始化时将参数设为mu,sigma
    :   __mu(mu),
        __sigma(sigma)
    {}
    double operator()(rand_16807& rand1, rand_16807& rand2)       //用box muller方法获得正态分布
    {
        double x, r, theta;
        r = sqrt(-2 * log(rand1() * 1.0 / rand_range));
        theta = 2.0 * M_PI * rand2() / rand_range;
        x = r * cos(theta);
        return x * __sigma + __mu;
    }
    void set_parameter(double mu, double sigma)       //将参数设为mu,sigma
    {
        __mu = mu;
        __sigma = sigma;
    }
private:
    double __gauss(double x)
    {
        return exp(-x * x / 2) / sqrt(2 * M_PI);
    }
    double __mu;
    double __sigma;
};

#endif // RANDOM_H_INCLUDED
