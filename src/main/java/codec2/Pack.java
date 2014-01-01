/*
 * Copyright (C) 2010 Perens LLC <bruce@perens.com>
 *
 * All Rights Reserved
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 2.1, as published
 * by the Free Software Foundation. This program is distributed in the hope that
 * it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 */
package codec2;

/**
 * Pack a long bit field into bytes, after first encoding the field with Gray code.
 *
 * <p>The output is an array of byte data, and the Gray coding is to reduce
 * the effect of single-bit errors. Although these bytes are not transmitted, they
 * are only exchanged with an application.  They could be transmitted in a higher
 * bandwidth modem.
 *
 * <p>Indices are always expected to be >= 0.
 *
 * <p>Copyright (C) 2010 Perens LLC<br>
 * All Rights Reserved
 */
public final class Pack {

    private static int bitOffset;

    public void Init() {
        bitOffset = 0;
    }

    public void pack(byte[] bitArray, int value, int valueBits) {
        /*
         * Convert from binary to Gray code
         */

        value = (value >> 1) ^ value;

        do {
            int bI = bitOffset;
            int bitsLeft = 8 - (bI & 0x7);
            int sliceWidth = bitsLeft < valueBits ? bitsLeft : valueBits;
            int wordIndex = bI >>> 3;

            bitArray[wordIndex] |= ((byte) ((value >> (valueBits - sliceWidth)) << (bitsLeft - sliceWidth)));

            bitOffset = bI + sliceWidth;
            valueBits -= sliceWidth;
        } while (valueBits != 0);
    }

    /*
     * Unpack a field from a byte array, converting from Gray code to binary.
     */
    public int unpack(byte[] bitArray, int valueBits) {
        int field = 0;
        int t;

        do {
            int bI = bitOffset;
            int bitsLeft = 8 - (bI & 0x7);
            int sliceWidth = bitsLeft < valueBits ? bitsLeft : valueBits;

            field |= (((bitArray[bI >>> 3] >> (bitsLeft - sliceWidth)) & ((1 << sliceWidth) - 1)) << (valueBits - sliceWidth));

            bitOffset = bI + sliceWidth;
            valueBits -= sliceWidth;
        } while (valueBits != 0);

        /*
         * Convert from Gray code to binary
         */

        t = field ^ (field >>> 8);
        t ^= (t >>> 4);
        t ^= (t >>> 2);
        t ^= (t >>> 1);

        return t;
    }
}
