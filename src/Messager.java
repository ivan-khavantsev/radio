import audio.AudioChannel;
import crypt.AES;
import modem.Demodulator;
import modem.Modulator;

import javax.sound.sampled.AudioInputStream;

public class Messager {
    private AES aes = null;

    public void start() {

        try {
            aes = new AES("abc");
        } catch (Throwable t) {
            System.out.println(t);
        }

        new Thread(new Runnable() {
            @Override
            public void run() {
                Demodulator demodulator = new Demodulator(new AudioInputStream(AudioChannel.getInLine()));
//                try {
//                   Thread.sleep(1000);
//                } catch(InterruptedException ex) {
//                    Thread.currentThread().interrupt();
//                }
                for (int i =0 ;i<10000;i++) {
                    try {
                        byte[] encrypted = demodulator.demodulate(48);
                        byte[] data = aes.decrypt(encrypted);
                        System.out.println(new String(data));
                    } catch (Throwable t) {
                        System.out.println(t);
                    }
                }
            }
        }).start();


        try {
            Modulator modulator = new Modulator();



            int i = 1000;
            while (true) {
                byte[] text = ("Hello my friend! How are you? . "+i).getBytes("UTF-8");

                byte[] encrypted = aes.encrypt(text);
                byte[] dataPacket = modulator.modulate(encrypted);
                AudioChannel.getOutLine().write(dataPacket, 0, dataPacket.length);
                i++;
            }

        } catch (Throwable t) {
            System.out.println(t);
        }

    }
}
