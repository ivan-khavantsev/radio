public class Main {


    public static void main(String[] args) throws Throwable {
        new Thread(new Sender()).start();
        new Thread(new Receiver()).start();
     //  new Messager().start();
    }


}
