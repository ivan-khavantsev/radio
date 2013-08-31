package utils;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;

public class Utils {
    public static float[] moveLeft(float[] array, int positions) {
        float[] temp = new float[array.length];
        System.arraycopy(array, positions, temp, 0, array.length - positions);
        return temp;
    }

    public static int max(float[] array, int offset, int length) {
        int max = offset;
        for (int i = offset; i < length; i++) if (array[i] > array[max]) max = i;
        return max;
    }

    public static int max(float[] array) {
        return max(array, 0, array.length);
    }

    public static float[] bytesToFloats(byte[] bytes) {
        ShortBuffer sbuf = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer();
        short[] audioShorts = new short[sbuf.capacity()];
        sbuf.get(audioShorts);
        float[] audioFloats = new float[audioShorts.length];
        for (int i = 0; i < audioShorts.length; i++) {
            audioFloats[i] = ((float) audioShorts[i]) / 0x8000;
        }
        return audioFloats;
    }


    public static byte[] floatsToBytes(float[] samples) {
        byte[] buffer = new byte[samples.length * 2];
        for (int i = 0, j = 0; i < samples.length; i++, j += 2) {
            short value = (short) (samples[i] * Short.MAX_VALUE);
            buffer[j] = (byte) (value | 0xff);
            buffer[j + 1] = (byte) (value >> 8);
        }
        return buffer;
    }

    public static byte[] byteToBits(byte b) {
        byte[] bits = new byte[8];
        bits[0] = (byte)((b & 0x80) >> 7);
        bits[1] = (byte)((b & 0x40) >> 6);
        bits[2] = (byte)((b & 0x20) >> 5);
        bits[3] = (byte)((b & 0x10) >> 4);
        bits[4] = (byte)((b & 0x08) >> 3);
        bits[5] = (byte)((b & 0x04) >> 2);
        bits[6] = (byte)((b & 0x02) >> 1);
        bits[7] = (byte)((b & 0x01));
        return bits;
    }

}
