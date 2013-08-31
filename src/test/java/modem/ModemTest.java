package modem;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import java.io.ByteArrayInputStream;

import static org.junit.Assert.assertEquals;

@RunWith(JUnit4.class)
public class ModemTest {
    @Test
    public void modem() throws Throwable {
        String text = "This is modem test";
        byte[] data = text.getBytes("UTF-8");

        Modulator modulator = new Modulator();

        byte[] modulatedData = modulator.modulate(data);
        ByteArrayInputStream inputStream = new ByteArrayInputStream(modulatedData);

        Demodulator demodulator = new Demodulator(inputStream);
        byte[] demodulatedData = demodulator.demodulate(data.length);
        String demodulatedText = new String(demodulatedData);

        assertEquals(text, demodulatedText);
    }
}
