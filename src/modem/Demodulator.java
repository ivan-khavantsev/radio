package modem;

import audio.analysis.FFT;
import audio.analysis.FourierTransform;
import utils.Utils;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

public class Demodulator {

    private InputStream inputStream;
    private boolean sync = false;

    public Demodulator(InputStream inputStream) {
        this.inputStream = new BufferedInputStream(inputStream, 320);
    }

    public synchronized byte[] demodulate(int packetSize) throws Exception {
        byte[] withoutSyncBytes = new byte[Modulator.SYNC_READY_BYTES.length];
        if (!sync) {
            this.sync();
            sync = true;
        } else {
            inputStream.read(withoutSyncBytes);
        }

        byte[] data = this.readBytes(packetSize);
        return data;
    }

    public void resync() {
        sync = false;
    }

    private void sync() throws IOException {
        byte[] syncBytes = new byte[Modulator.SYNC_READY_BYTES.length];
        inputStream.read(syncBytes, 0, syncBytes.length);
        float[] audioFloats = Utils.bytesToFloats(syncBytes);

        while (true) {
            float syncCorelation = getCorrelation(audioFloats, Modulator.SYNC_SAMPLES);
            if (syncCorelation > 0.8f) {
                break;
            } else {
                audioFloats = Utils.moveLeft(audioFloats, 1);
                byte[] tempSampleBytes = new byte[Modulator.BYTES_PER_AUDIO_SAMPLE];
                inputStream.read(tempSampleBytes, 0, tempSampleBytes.length);
                float[] tempSample = Utils.bytesToFloats(tempSampleBytes);
                audioFloats[audioFloats.length - 1] = tempSample[0];
            }
        }

    }


    private byte[] readBytes(int count) throws IOException {
        byte[] result = new byte[count];
        byte[] buffer = new byte[Modulator.SAMPLES_PER_DATA_BIT * Modulator.BYTES_PER_AUDIO_SAMPLE];
        float[] audioFloats;
        for (int i = 0; i < count; i++) {
            byte b = 0;
            for (int j = 8; j > 0; j--) {
                inputStream.read(buffer, 0, buffer.length);
                audioFloats = Utils.bytesToFloats(buffer);
                int bit = this.samplesToBit(audioFloats);
                b = (byte) (b | bit << (j - 1));
            }
            result[i] = b;
        }
        return result;
    }

    private int samplesToBit(float[] samples) {
        float correlation = getCorrelation(samples, Modulator.ZERO_SAMPLES);
        if (correlation > 0) {
            return 0;
        } else {
            return 1;
        }
    }

    //TODO: CALC FOR POPULAR ETALON
    /*
     * Корреляция Пирсона
     * http://cito-web.yspu.org/link1/metod/met125/node35.html
     */
    public static float getCorrelation(float[] signal1, float[] signal2) {
        float signal1Sum = sum(signal1, false);
        float signal2Sum = sum(signal2, false);

        float signalsMultSum = 0;
        for (int i = 0; i < signal1.length; i++) {
            signalsMultSum += signal1[i] * signal2[i];
        }
        float top = (signal1.length * signalsMultSum) - (signal1Sum * signal2Sum);
        float signal1QuadSum = sum(signal1, true);
        float signal2QuadSum = sum(signal2, true);

        float bottom1 = signal1.length * signal1QuadSum - signal1Sum * signal1Sum;
        float bottom2 = signal2.length * signal2QuadSum - signal2Sum * signal2Sum;

        float bottom = (float) Math.sqrt(bottom1 * bottom2);

        float result = top / bottom;
        return result;
    }

    public static float sum(float[] values, boolean quad) {
        float sum = 0;
        for (float v : values) {
            if (quad) {
                sum += v * v;
            } else {
                sum += v;
            }
        }
        return sum;
    }

}
