import os
import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QLabel, QPushButton, QLineEdit, QTextEdit, QVBoxLayout, QWidget, QMessageBox, QDialog, QDialogButtonBox, QVBoxLayout, QFormLayout
from datetime import datetime
import requests
import getpass
from cryptography.fernet import Fernet, InvalidToken

class CredentialsDialog(QDialog):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("LASP WebIAM Credentials")

        layout = QVBoxLayout()

        self.username_entry = QLineEdit()
        self.password_entry = QLineEdit()
        self.password_entry.setEchoMode(QLineEdit.EchoMode.Password)

        form_layout = QFormLayout()
        form_layout.addRow("LASP WebIAM Username:", self.username_entry)
        form_layout.addRow("LASP WebIAM Password:", self.password_entry)

        layout.addLayout(form_layout)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        layout.addWidget(button_box)

        self.setLayout(layout)

    def get_credentials(self):
        username = self.username_entry.text()
        password = self.password_entry.text()

        return username, password

class DataDownloaderApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Data Downloader")
        self.setGeometry(100, 100, 400, 300)
        self.setup_ui()

        # Initialize key attribute
        self.key = self.load_key()
        # Load credentials if available
        self.load_credentials()

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout()

        start_date_label = QLabel("Start Date (YYYY-MM-DD):")
        layout.addWidget(start_date_label)
        self.start_date_entry = QLineEdit("2024-05-09")
        layout.addWidget(self.start_date_entry)

        stop_date_label = QLabel("Stop Date (YYYY-MM-DD):")
        layout.addWidget(stop_date_label)
        self.stop_date_entry = QLineEdit("2024-05-10")
        layout.addWidget(self.stop_date_entry)

        self.download_button = QPushButton("Download Data (.bin Format)")
        self.download_button.clicked.connect(self.download_data)
        layout.addWidget(self.download_button)

        self.get_info_button = QPushButton("Get Info")
        self.get_info_button.clicked.connect(self.get_info)
        self.get_info_button.setEnabled(True)
        layout.addWidget(self.get_info_button)

        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        layout.addWidget(self.info_text)

        central_widget.setLayout(layout)

        # Initialize URL parameters
        self.url_base = "https://lasp.colorado.edu/ops/imap/poda/dap2/packets/SID2/IDX_SCI"
        self.url_project = "project(packet)"

    def download_data(self):
        start_date = self.start_date_entry.text()
        stop_date = self.stop_date_entry.text()

        try:
            datetime.strptime(start_date, "%Y-%m-%d")
            datetime.strptime(stop_date, "%Y-%m-%d")
        except ValueError:
            QMessageBox.critical(self, "Error", "Invalid date format. Please use YYYY-MM-DD.")
            return

        # Construct the URL
        url = f"{self.url_base}.bin?time>={start_date}T00:00:01&time<{stop_date}T18:27:00&{self.url_project}&formatTime(%22yyyy-MM-dd%27T%27HH:mm:ss%22)"

        # Add authentication headers
        username = self.username
        password = self.password
        auth = requests.auth.HTTPBasicAuth(username, password)

        # Download the data
        response = requests.get(url, auth=auth)

        if response.status_code == 200:
            # Write the binary data to a local file
            file_path = f"IDX_SCI_{start_date}_to_{stop_date}.bin"
            with open(file_path, "wb") as f:
                f.write(response.content)
            QMessageBox.information(self, "Success", f"Data downloaded successfully and saved to '{file_path}'")
            self.get_info_button.setEnabled(True)
        else:
            QMessageBox.critical(self, "Error", f"Failed to download data. Status code: {response.status_code}")
            self.get_info_button.setEnabled(False)

    def load_credentials(self):
        if os.path.exists("credentials.enc"):
            with open("credentials.enc", "rb") as file:
                encrypted_data = file.read()
                print("||===|| Credentials read ||===||", encrypted_data)
            fernet = Fernet(self.key)
            try:
                decrypted_data = fernet.decrypt(encrypted_data).decode().splitlines()
                print("Decrypted data:", decrypted_data)
                if len(decrypted_data) == 2:
                    self.username, self.password = decrypted_data
            except InvalidToken:
                print("Invalid token or key. Credentials decryption failed.")
        else:
            self.get_credentials()


    def get_credentials(self):
        dialog = CredentialsDialog()
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.username, self.password = dialog.get_credentials()
            save_credentials = QMessageBox.question(self, "Save Credentials", "Do you want to save these credentials?", QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if save_credentials == QMessageBox.StandardButton.Yes:
                self.save_credentials(self.username, self.password)
        else:
            QMessageBox.critical(self, "Error", "Credentials input canceled.")

    def save_credentials(self, username, password):
        if not self.key:
            self.key = Fernet.generate_key()
            self.save_key()

        fernet = Fernet(self.key)
        encrypted_data = fernet.encrypt(f"{username}\n{password}".encode())

        with open("credentials.enc", "wb") as file:
            file.write(encrypted_data)

    def get_info(self):
        start_date = self.start_date_entry.text()
        stop_date = self.stop_date_entry.text()

        self.info_text.clear()

        # Display information in the text box
        self.info_text.append(f"Start Date: {start_date}")
        self.info_text.append(f"Stop Date: {stop_date}")

        # Get URL for .ascii format
        url_ascii = f"https://lasp.colorado.edu/ops/imap/poda/dap2/packets/SID2/IDX_SCI.asc?time%3E={start_date}T00:00:01&time%3C{stop_date}T18:27:00&project(time,packet)&formatTime(%22yyyy-MM-dd%27T%27HH:mm:ss%22)"

        # Add authentication credentials
        username = self.username
        password = self.password
        
        # Download the data
        response = requests.get(url_ascii, auth=(username, password))

        if response.status_code == 200:
            # Write the ASCII data to a local file
            file_path = f"IDX_SCI_{start_date}_to_{stop_date}.asc"
            with open(file_path, "w") as f:
                f.write(response.text)
            total_bytes = 0  # Initialize total bytes counter
            packetnum = 0

            # Read the contents of the file and print to self.info_text
            with open(file_path, "r") as f:
                ascii_content = f.read()
                self.info_text.append("Downloaded Data (.ascii Format):")
                self.info_text.append(ascii_content)
                
                # Reset file pointer to the beginning of the file
                f.seek(0)
                
                for line in f:
                    if "Byte Binary" in line:
                        packetnum += 1
                        byte_amount = int(line.split()[2])  # Extract byte amount from the line
                        total_bytes += byte_amount  # Update total bytes counter
                    self.info_text.append(line.strip())  # Add each line to the text box

                # Display number of packets transmitted              
                self.info_text.append(f"Science Packets Transmitted: {packetnum} packets")
                # Display total bytes transmitted
                self.info_text.append(f"Total Bytes Transmitted: {total_bytes} bytes")
                    
            QMessageBox.information(self, "Success", f"Data downloaded successfully and saved to '{file_path}'")
        else:
            QMessageBox.critical(self, "Error", f"Failed to download data. Status code: {response.status_code}")

    def save_key(self):
        with open("encryption_key.txt", "wb") as key_file:
            key_file.write(self.key)

    def load_key(self):
        if os.path.exists("encryption_key.txt"):
            with open("encryption_key.txt", "rb") as key_file:
                return key_file.read()
        else:
            return None


def main():
    app = QApplication(sys.argv)
    QApplication.setStyle("fusion")  # Set application style to "fusion"

    window = DataDownloaderApp()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
