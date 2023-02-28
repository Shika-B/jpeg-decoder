#include "utils.h"

unsigned short u8_to_u16(unsigned char x, unsigned char y) {
    return (x << 8) + y;
}