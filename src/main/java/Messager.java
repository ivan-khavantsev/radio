import audio.AudioChannel;
import crypt.AES;
import modem.Demodulator;
import modem.Modulator;
import reedsolomon.GenericGF;
import reedsolomon.ReedSolomonDecoder;
import reedsolomon.ReedSolomonEncoder;

import javax.sound.sampled.AudioInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;

public class Messager {
    private AES aes = null;

    public void start() {
        try {
            aes = new AES("abc");
        } catch (Throwable t) {
            System.out.println(t);
        }

        GenericGF gf = GenericGF.DATA_MATRIX_FIELD_256;


        ReedSolomonEncoder solomonEncoder = new ReedSolomonEncoder(gf);
        final ReedSolomonDecoder solomonDecoder = new ReedSolomonDecoder(gf);


        new Thread(new Runnable() {
            @Override
            public void run() {
                Demodulator demodulator = new Demodulator(new AudioInputStream(AudioChannel.getInLine()));
//                try {
//                   Thread.sleep(1000);
//                } catch(InterruptedException ex) {
//                    Thread.currentThread().interrupt();
//                }
                for (int i = 0; i < 10000; i++) {
                    try {
                        byte[] encrypted = demodulator.demodulate(18);
                        int[] solomoned = new int[18];
                        ByteToInt(solomoned,encrypted,18);
                        solomonDecoder.decode(solomoned, 2);
                        byte[] desolomoned = new byte[16];
                        IntToByte(desolomoned,solomoned,16);



                        byte[] data = aes.decrypt(desolomoned);
                        System.out.println(new String(data));
                    } catch (Throwable t) {
                        System.out.println("err1: " + t);
                        demodulator.resync();
                    }
                }
            }
        }).start();


        Modulator modulator = new Modulator();
        int i = 1000;
        while (true) {
            try {
                byte[] text = (i + " ").getBytes("UTF-8");
                byte[] encrypted = aes.encrypt(text);

                int[] solomoned = new int[18];
                ByteToInt(solomoned,encrypted,16);
                solomonEncoder.encode(solomoned, 2);



                byte[] solomonedBytes = new byte[18];
                IntToByte(solomonedBytes, solomoned, 18);


                byte[] dataPacket = modulator.modulate(solomonedBytes);
                AudioChannel.getOutLine().write(dataPacket, 0, dataPacket.length);
                i++;
            } catch (Throwable t) {
                System.out.println("err2: " + t);
            }
        }

    }


    public void IntToByte(byte arrayDst[], int arrayOrg[], int maxOrg){
        for (int i=0; i<maxOrg;i++){
            arrayDst[i] = (byte)arrayOrg[i];
        }
    }

    public void ByteToInt(int arrayDst[], byte arrayOrg[], int maxOrg){
        for (int i=0; i<maxOrg;i++){
            arrayDst[i] = arrayOrg[i] & 0xff;
        }
    }



}
