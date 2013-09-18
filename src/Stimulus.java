import javax.swing.JOptionPane;


public class Stimulus {
double duration,applied_t,voltage;
    Stimulus(){
        duration=0.0;
        applied_t=0.0;
        voltage=0.0;
    }
    public void set(){
         String userInput = JOptionPane.showInputDialog("Enter Duration");
         duration = Double.parseDouble(userInput);
         userInput = JOptionPane.showInputDialog("Enter time to apply Stimulus");
         applied_t = Double.parseDouble(userInput);
         userInput = JOptionPane.showInputDialog("Enter Voltage");
         voltage = Double.parseDouble(userInput);
    }
    public static void main(String[] args) {
        Stimulus b= new Stimulus();
        b.set();
    }
}
