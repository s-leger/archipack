
import array

"""
    basic bitarray
"""
class BitArray():
    
    def __init__(self, bitSize, fill = 0):
        self.size = bitSize
        intSize = bitSize >> 5
        if (bitSize & 31):
            intSize += 1
        if fill == 1:
            fill = 4294967295
        else:
            fill = 0
        self.bitArray = array.array('I')
        self.bitArray.extend((fill,) * intSize)

    def __str__(self):
        return str(self.list)

    def bit_location(self, bit_num):
        return bit_num >> 5, bit_num & 31

    def test(self, bit_num):
        record, offset = self.bit_location(bit_num)
        mask = 1 << offset
        return(self.bitArray[record] & mask)

    def set(self, bit_num):
        record, offset = self.bit_location(bit_num)
        mask = 1 << offset
        self.bitArray[record] |= mask

    def clear(self, bit_num):
        record, offset = self.bit_location(bit_num)
        mask = ~(1 << offset)
        self.bitArray[record] &= mask

    def toggle(self, bit_num):
        record, offset = self.bit_location(bit_num)
        mask = 1 << offset
        self.bitArray[record] ^= mask

    @property 
    def len(self):
        return len(self.bitArray)
        
    @property
    def copy(self):
        copy = BitArray(self.size)
        for i in range(self.len):
            copy.bitArray[i] = self.bitArray[i]
        return copy        
    
    @property
    def list(self):
        return [x for x in range(self.size) if self.test(x) > 0]

    def none(self):
        for i in range(self.len):
            self.bitArray[i] = 0

    def reverse(self):
        for i in range(self.len):
            self.bitArray[i] = 4294967295 ^ self.bitArray[i]

    def all(self):
        for i in range(self.len):
            self.bitArray[i] = 4294967295
    

        