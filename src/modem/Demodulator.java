package modem;

import audio.analysis.FFT;
import audio.analysis.FourierTransform;
import utils.Plot;

import java.awt.*;
import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;
import java.util.LinkedList;
import java.util.Queue;

public class Demodulator {

    class InBuffer extends InputStream {
        private InputStream is = null;
        private ByteBuffer bb = null;
        Queue<Byte> arrayQueue = new LinkedList<Byte>();

        public InBuffer(InputStream is) {
            this.is = new BufferedInputStream(is, 2);
        }


        public void add(byte[] data, int offset, int count) {
            for (int i = offset; i < offset + count; i++) {
                arrayQueue.add(data[i]);
            }
        }


        @Override
        public int read() throws IOException {
            int ret = -1;
            if (!arrayQueue.isEmpty()) {
                ret = arrayQueue.poll() & 0xFF;
            } else {
                ret = is.read();
            }
            return ret;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
    //Plot plot = new Plot("Note A Spectrum", 1280, 512);
    private InBuffer inputStream;
    private FourierTransform fft = new FFT(Modulator.SAMPLES_PER_BIT, Modulator.SAMPLE_RATE);
    // private FourierTransform syncFft = new FFT(SYNC_SAMPLES_PER_SYMBOL, Modulator.SAMPLE_RATE);

    public Demodulator(InputStream inputStream) {
        this.inputStream = new InBuffer(inputStream);
    }

    /*
        @int packetSize - size of data packet to demodulate
     */
    public synchronized byte[] demodulate(int packetSize) throws Exception {
        ByteBuffer dataBuffer = ByteBuffer.allocate(packetSize);
        this.sync1(packetSize);
        for (int i = 0; i < packetSize; i++) {
            byte b = this.readByte();
            dataBuffer.put(b);
        }
        return dataBuffer.array();
    }

//    private void sync() throws IOException {
//        int[] syncTrue = new int[4];
//        while (true) {
//            byte[] syncBytes = new byte[SYNC_SAMPLES_PER_SYMBOL * 2];
//            inputStream.read(syncBytes, 0, syncBytes.length);
//            float[] syncFloats = getFloatSamples(syncBytes);
//            syncFft.forward(syncFloats);
//            int max = max(syncFft.getSpectrum(), 1, 4);
//            syncTrue[syncTrue.length-1] = max;
//            if(syncTrue[0]==2 && syncTrue[1]==1 && syncTrue[2]==4 && syncTrue[3] == 3){
//                break;
//            }else{
//                syncTrue = moveLeft(syncTrue,1);
//            }
//        }
//
//
//    }

    private void sync1(int packetSize) throws IOException {
        /*
            Синхронизация последовательностью Бэкера
            sync = [+1, +1, +1, -1, -1, +1, -1];
         */
        //((packetSize * Modulator.SAMPLES_PER_BYTE)*Modulator.BYTES_PER_SAMPLE) +
        int fullPacketSize = (Modulator.SYNC_SEQUENCE.length * Modulator.SYNC_SAMPLES_PER_SYMBOL) * 2;
        byte[] syncBytes = new byte[fullPacketSize];
        inputStream.read(syncBytes, 0, syncBytes.length);
        float[] audioFloats = this.getFloatSamples(syncBytes);


        while (true) {
            int ii = max(audioFloats);
            float max = audioFloats[ii];



            if (max > 0.3f && rf(audioFloats) > max-0.1f) {
                break;
            } else {
                audioFloats = moveLeft(audioFloats, 1);
                byte[] tmp = new byte[2];
                inputStream.read(tmp, 0, tmp.length);
                float[] dochitannieSamples = getFloatSamples(tmp);
                audioFloats[audioFloats.length - 1] = dochitannieSamples[0];
            }
        }
    }



    public float[] getRF(float[] samples) {
        float[] rf = new float[samples.length]; // ??? - Значения взаимной корреляционной функции
        float[] sums = new float[Modulator.SYNC_SEQUENCE.length];
        float r = 0;
        for (int i = 0; i < Modulator.SYNC_SEQUENCE.length - 1; i++) {
            float s = 0;
            for (int j = 0; j < Modulator.SYNC_SEQUENCE.length - 1; j++) {
                s = s + samples[i * Modulator.SYNC_SAMPLES_PER_SYMBOL + j];
            }
            sums[i] = s;
            r = r + (sums[i] * Modulator.SYNC_SEQUENCE[i]) / Modulator.SYNC_SAMPLES_PER_SYMBOL;
        }

        rf[0] = r / Modulator.SYNC_SEQUENCE.length;


        for (int i = 2; i < (samples.length - Modulator.SYNC_SAMPLES_PER_SYMBOL * Modulator.SYNC_SEQUENCE.length); i++) {
            r = 0;
            for (int j = 1; j < Modulator.SYNC_SEQUENCE.length; j++) {
                sums[j] = sums[j] - samples[i - 1 + (j - 1) * Modulator.SYNC_SAMPLES_PER_SYMBOL] + samples[i - 1 + j * Modulator.SYNC_SAMPLES_PER_SYMBOL];
                r = r + sums[j] * Modulator.SYNC_SEQUENCE[j] / Modulator.SYNC_SAMPLES_PER_SYMBOL;
            }
            rf[i] = r / Modulator.SYNC_SEQUENCE.length;
        }
        //plot.clear();
        // plot.plot(samples, 1, Color.blue);
        //  plot.plot(rf, 1, Color.red);
        return rf;
    }

    public static float[] moveLeft(float[] array, int positions) {
        float[] temp = new float[array.length];
        System.arraycopy(array, positions, temp, 0, array.length - positions);
        return temp;
    }

    public static int[] moveLeft(int[] array, int positions) {
        int[] temp = new int[array.length];
        System.arraycopy(array, positions, temp, 0, array.length - positions);
        return temp;
    }


    private byte readByte() throws IOException {
        int leftForByte = 8;
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < 8; i++) {
            byte[] buffer = new byte[Modulator.SAMPLES_PER_BIT * Modulator.BYTES_PER_SAMPLE];
            inputStream.read(buffer, 0, buffer.length);
            float[] audioFloats = this.getFloatSamples(buffer);
            int bit = this.getSamplesBit(audioFloats);

            sb.append(bit);

            leftForByte--;

            if (leftForByte < 1) {
                byte[] bval = new BigInteger(sb.toString(), 2).toByteArray();
                if (bval.length == 2) {
                    return bval[1];
                } else {
                    return bval[0];
                }
            }
        }
        return 0;
    }

    private synchronized int getSamplesBit(float[] audioFloats) {
        fft.forward(audioFloats);
        float[] spectrum = fft.getSpectrum();

        int bit;
        if (spectrum[0] > spectrum[1]) {
            bit = 0;
        } else {
            bit = 1;
        }
        return bit;
    }


    public static float[] getFloatSamples(byte[] bytes) {
        ShortBuffer sbuf = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer();
        short[] audioShorts = new short[sbuf.capacity()];
        sbuf.get(audioShorts);
        float[] audioFloats = new float[audioShorts.length];
        for (int i = 0; i < audioShorts.length; i++) {
            audioFloats[i] = ((float) audioShorts[i]) / 0x8000;
        }
        return audioFloats;
    }

    public static int max(float[] array, int offset, int length) {
        int max = offset;
        for (int i = offset; i < length; i++) if (array[i] > array[max]) max = i;
        return max;
    }

    public static int max(float[] array) {
        int max = 0;
        for (int i = 0; i < array.length; i++) if (array[i] > array[max]) max = i;
        return max;
    }


    public float rf(float[] samples) {
        float[] sums = new float[Modulator.SYNC_SEQUENCE.length];
        float r = 0;
        for (int i = 0; i < Modulator.SYNC_SEQUENCE.length; i++) {
            float s = 0;
            for (int j = 0; j < Modulator.SYNC_SAMPLES_PER_SYMBOL; j++) {
                int cc = i * Modulator.SYNC_SAMPLES_PER_SYMBOL + j;
                s = s + samples[cc];
            }
            sums[i] = s;
            r = r + (sums[i] * Modulator.SYNC_SEQUENCE[i]) / Modulator.SYNC_SAMPLES_PER_SYMBOL;
        }

        float result = r / Modulator.SYNC_SEQUENCE.length;
        return result;
    }

}
