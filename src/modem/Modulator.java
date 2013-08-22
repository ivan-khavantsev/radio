package modem;

import utils.Utils;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;

public class Modulator {

    public static final int SAMPLE_RATE = 44100;
    public static final int SAMPLES_PER_DATA_BIT = 16;
    public static final int SAMPLES_PER_DATA_BYTE = SAMPLES_PER_DATA_BIT * 8;
    public static final int BYTES_PER_AUDIO_SAMPLE = 2;

    public static final byte[] SYNC_READY_BYTES;
    public static final float[] SYNC_SEQUENCE = {+1, +1, +1, -1, -1, +1, -1};
    public static final int SYNC_SAMPLES_PER_SYMBOL = 4;

    private static final float[] FREQUENCIES = {689, 2756};   //0 - 689 MHz, 1 -2756 MHz
    private float angle = 0;
    public static final float[] SYNC_SAMPLES;

    static {
//        float[] tempSync = new float[SYNC_SEQUENCE.length * SYNC_SAMPLES_PER_SYMBOL];
//        for (int i = 0; i < SYNC_SEQUENCE.length; i++) {
//            for (int j = 0; j < SYNC_SAMPLES_PER_SYMBOL; j++) {
//                tempSync[i * SYNC_SAMPLES_PER_SYMBOL + j] = SYNC_SEQUENCE[i];
//            }
//        }

        float[] tempSync = new float[SAMPLES_PER_DATA_BIT*2];

        float[] upHalfWave = getWave(8,0);
        float[] downHalfWave = getWave(8,Math.PI);

        System.arraycopy(upHalfWave, 0, tempSync, 0, upHalfWave.length);
        System.arraycopy(upHalfWave, 0, tempSync, 8, upHalfWave.length);
        System.arraycopy(downHalfWave, 0, tempSync, 16, downHalfWave.length);
        System.arraycopy(downHalfWave, 0, tempSync, 24, downHalfWave.length);
        SYNC_SAMPLES = tempSync;
        SYNC_READY_BYTES = Utils.floatsToBytes(tempSync);
    }


    public static final float[] ZERO_SAMPLES;
    public static final float[] ONE_SAMPLES;

    static {
       ZERO_SAMPLES = getWave(SAMPLES_PER_DATA_BIT, 0);
       ONE_SAMPLES = getWave(SAMPLES_PER_DATA_BIT, Math.PI);
    }


    public byte[] modulate(byte[] data) {
        ByteBuffer samplesByteBuffer = ByteBuffer.allocate(
                (data.length * SAMPLES_PER_DATA_BYTE * BYTES_PER_AUDIO_SAMPLE) + SYNC_READY_BYTES.length);
        samplesByteBuffer.put(SYNC_READY_BYTES);

        for (byte b : data) {
            float[] floatSamples = this.getByteSamples(b);
            byte[] byteSamples = Utils.floatsToBytes(floatSamples);
            samplesByteBuffer.put(byteSamples);
        }

        angle = 0;
        return samplesByteBuffer.array();
    }

    private float[] getByteSamples(byte b) {
        FloatBuffer floatBuffer = FloatBuffer.allocate(SAMPLES_PER_DATA_BYTE);
        byte[] bits = Utils.byteToBits(b);
        for (byte bit : bits) {
            floatBuffer.put(getBitSamples(bit));
        }
        float[] buffer = floatBuffer.array();
        return buffer;
    }

    private float[] getBitSamples1(byte bit) {
        float frequency = FREQUENCIES[bit];
        float increment = (float) (2 * Math.PI) * frequency / SAMPLE_RATE;
        float samples[] = new float[SAMPLES_PER_DATA_BIT];
        for (int i = 0; i < samples.length; i++) {
            samples[i] = (float) Math.sin(angle);
            angle += increment;
        }
        return samples;
    }

    private float[] getBitSamples(byte bit) {
       if(bit == 0){
           return ZERO_SAMPLES;
       }else {
           return ONE_SAMPLES;
       }
    }




    public static final float FREQUENCY1 = 2756.25f;
    public static float[] getWave(int samples, double phaseAngle) {
        float[] wave = new float[samples];
        double angle1 = phaseAngle;
        float increment1 = (float) (2 * Math.PI) * FREQUENCY1 / 44100;
        for (int i = 0; i < wave.length; i++) {
            wave[i] =  (float)( Math.sin(angle1));
           // if(noise) wave[i] += Math.random()/2;
            angle1 += increment1;
        }
        return wave;
    }

}
