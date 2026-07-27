#include "ui/MainWindow.h"

#include <QApplication>
#include <QFont>

int main(int argc, char* argv[]) {
    QApplication application(argc, argv);
    application.setApplicationName(QStringLiteral("Stochia"));
    application.setApplicationDisplayName(QStringLiteral("Stochia · 随机空间"));
    application.setOrganizationName(QStringLiteral("Stochia Lab"));
    application.setApplicationVersion(QStringLiteral("0.2.1"));
    application.setStyle(QStringLiteral("Fusion"));
    application.setFont(QFont(QStringLiteral("Segoe UI"), 10));

    MainWindow window;
    window.show();
    return application.exec();
}
