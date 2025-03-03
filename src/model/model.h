#ifndef MODEL_H
#define MODEL_H


// probability precision
#define PBITS   (15)
#define PSCALE  (1<<PBITS)
#define PSCALEh (PSCALE>>1)
#define PSCALEm (PSCALE-1)

// weight precision
#define WBITS   (16)
#define WSCALE  (1<<WBITS)
#define WSCALEh (WSCALE>>1)

template <class T> T clamp(T val,T min,T max) {return val<min?min:val>max?max:val;};

#endif // MODEL_H
