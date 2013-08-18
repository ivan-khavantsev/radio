package audio.samples.part5;

import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;

import audio.analysis.FFT;
import org.xiph.speex.SpeexDecoder;

import javax.sound.sampled.*;

/**
 * Simple example that generates a 1024 samples sine wave at 440Hz
 * and plots the resulting spectrum.
 *
 * @author mzechner
 */
public class FourierTransformPlot {


    public static byte b;
    public static StringBuilder sb = new StringBuilder();


    public static float[] floatFromBytes(byte[] bytes) {
        ShortBuffer sbuf =
                ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer();
        short[] audioShorts = new short[sbuf.capacity()];
        sbuf.get(audioShorts);
        float[] audioFloats = new float[audioShorts.length];
        for (int i = 0; i < audioShorts.length; i++) {
            audioFloats[i] = ((float) audioShorts[i]) / 0x8000;
        }
        return audioFloats;
    }

    static SpeexDecoder speexDecoder = new SpeexDecoder();
    static SourceDataLine line2 = null;

    public static void main(String[] argv) throws Throwable {

        AudioFormat audioFormat = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED,
                44100.0F, 16, 1, 2, 44100.0F, false);

        final TargetDataLine line = (TargetDataLine) AudioSystem.getLine(new DataLine.Info(TargetDataLine.class, audioFormat));
        line.open();
        line.start();


        AudioFormat audioFormat2 = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED,
                8000.0F, 16, 1, 2, 8000.0F, false);
        Mixer.Info[] mixerInfo = AudioSystem.getMixerInfo();
        Mixer mixer = AudioSystem.getMixer(mixerInfo[2]);
        line2 = (SourceDataLine) mixer.getLine(new DataLine.Info(SourceDataLine.class, audioFormat2));
        line2.open();
        line2.start();


        speexDecoder.init(0, 8000, 1, true);


        byte[] buffer = new byte[2 * 16];


        sync(line);


        while (true) {
            line.read(buffer, 0, buffer.length);

            float[] audioFloats = floatFromBytes(buffer);

            FFT fft = new FFT(audioFloats.length, 44100);
            fft.forward(audioFloats);
            float[] spect = fft.getSpectrum();
            if ((spect[0] - spect[1]) > 2f || (spect[1] - spect[0]) > 2f) {
                int bit = 0;
                if (spect[0] > spect[1]) {
                    bit = 0;
                } else {
                    bit = 1;
                }
                handleBit(bit);
            } else {
                //синхронизация
                //System.out.println();
            }

        }


//        final float frequency = 1000; // Note A
//        float increment = (float) (2 * Math.PI) * frequency / 44100;
//        float angle = 0;
//        float samples[] = new float[16];


//        for (int i = 0; i < samples.length; i++) {
//            samples[i] = (float) Math.sin(angle);
//            angle += increment;
//        }


        // for (int i = 0; i < 33; i++) {
//            System.out.println(i + ": " + fft.indexToFreq(i));
//        }
//
//        utils.Plot plot = new utils.Plot("Note A Spectrum", 512, 512);
//        plot.plot(fft.getSpectrum(), 1, Color.red);
    }

    private static void sync(TargetDataLine line) {
        byte[] syncBytes = new byte[2];
        int zeros = 0;
        while (true) {
            line.read(syncBytes, 0, syncBytes.length);
            float[] audioFloats = floatFromBytes(syncBytes);
            System.out.println(audioFloats[0]);
            if (audioFloats[0] < 0.02f && audioFloats[0] > -0.02f) {
                zeros++;
            } else {
                if (zeros > 5 && audioFloats[0] > 0.07f) {
                    break;
                } else {
                    zeros = 0;
                }
            }
        }

    }


    public static int sdvig = 7;
    public static ByteBuffer encoded = ByteBuffer.allocate(10);

    public static void handleBit(int bit) throws Throwable {
        if (sdvig < 0) {
            byte[] bval = new BigInteger(sb.toString(), 2).toByteArray();
            byte[] data = null;
            if(bval.length == 2){
                data = new byte[]{bval[1]};
            } else{
                data = new byte[]{bval[0]};
            }
            encoded.put(data);
            if (encoded.remaining() == 1) {
                byte[] encodedArray = encoded.array();
                speexDecoder.processData(encodedArray, 0, encodedArray.length);
                byte[] decodedData = new byte[speexDecoder.getProcessedDataByteSize()];
                speexDecoder.getProcessedData(decodedData, 0);
                line2.write(decodedData, 0, decodedData.length);


                encoded = ByteBuffer.allocate(320);
            }


            //System.out.print(new String(bval));
            //System.out.print(sb.toString());
            sdvig = 7;
            sb = new StringBuilder();


        }
        sb.append(bit);
        sdvig--;
    }
}
